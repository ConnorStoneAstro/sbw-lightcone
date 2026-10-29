"""Gallery of DES DR2 coadd cutouts (NOIRLab Data Lab SIA, des_dr2) for des_disagreement_cases, drawn with the
panels of src/dp02_gallery.py: the SBW most probable host with its half-light ellipse (solid), the DLR host with
its DLR ellipse (dashed), both at --ellipse-scale x their axes, and the published DES-SN host position.

For each case the SIA service lists the DR2 tiles covering the SN (searching 0.5 deg when its footprints
return nothing, which happens at some tile edges); the band's image (``<tile>_<band>.fits.fz``,
extension 1) of the tile whose centre is nearest is cut out through its access_url with POS / SIZE set to the
requested square. Requests are throttled (--min-interval) and cached under --cache-dir."""
from __future__ import annotations

import argparse
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyvo
import requests

import dp02_gallery as G

SIA = "https://datalab.noirlab.edu/sia/des_dr2"
G.MARKER_STYLES["pub"] = ("D", "black", 17, "published DES-SN host")
G.TITLE = "DES DR2 coadd"


def cutout_url(ra, dec, band, size, sia):
    t = sia.search((ra, dec), size=size / 3600.0).to_table()
    if len(t) == 0:   # the SIA footprints miss some tile edges: widen the search, keep the nearest tile centre
        t = sia.search((ra, dec), size=0.5).to_table()
    url = np.array([str(u) for u in t["access_url"]])
    ok = np.array([bool(re.search(rf"_{band}\.fits\.fz&extn=1&", u)) for u in url])
    if not ok.any():
        raise RuntimeError(f"no DR2 {band} tile covers ({ra:.5f}, {dec:.5f})")
    d = np.hypot((np.asarray(t["s_ra"], float) - ra) * np.cos(np.deg2rad(dec)), np.asarray(t["s_dec"], float) - dec)
    k = np.flatnonzero(ok)[np.argmin(d[ok])]
    s = f"{size / 3600.0:.8f}"
    return re.sub(r"&POS=[^&]*&SIZE=[^&]*", f"&POS={ra:.7f},{dec:.7f}&SIZE={s},{s}", url[k])


def fetch(ra, dec, band, size, sia, session, cache_dir, limiter, retries=4):
    path = cache_dir / f"des_{band}_{ra:.6f}_{dec:+.6f}_{size:g}as.fits"
    if path.exists():
        return path.read_bytes()
    err = None
    for attempt in range(retries):
        limiter.wait()
        try:
            r = session.get(cutout_url(ra, dec, band, size, sia), timeout=120)
            r.raise_for_status()
            if not r.content.startswith(b"SIMPLE"):
                raise ValueError(f"not FITS: {r.content[:80]!r}")
            cache_dir.mkdir(parents=True, exist_ok=True)
            path.write_bytes(r.content)
            return r.content
        except Exception as e:  # noqa: BLE001 - retry anything transient
            err = e
            time.sleep(3.0 * 2 ** attempt)
    raise RuntimeError(f"cutout failed after {retries} tries: {err}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cases")
    p.add_argument("--output", required=True)
    p.add_argument("--band", default="r", choices=list("grizy"))
    p.add_argument("--size", type=float, default=30.0)
    p.add_argument("--ncol", type=int, default=4)
    p.add_argument("--per-page", type=int, default=12)
    p.add_argument("--cache-dir", default="des_cutouts")
    p.add_argument("--ellipse-scale", type=float, default=1.0)
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--min-interval", type=float, default=0.7)
    p.add_argument("--dpi", type=int, default=150)
    args = p.parse_args(argv)
    for k, v in dict(stretch="asinh", pmin=20.0, pmax=99.7, dark=False, scalebar=5.0).items():
        setattr(args, k, v)
    cases = pd.read_csv(args.cases)
    sia, session = pyvo.dal.SIAService(SIA), requests.Session()
    limiter, cache = G.RateLimiter(args.min_interval), Path(args.cache_dir)
    images, errors = {}, {}

    def job(case):
        return G.read_image(fetch(float(case.sn_ra), float(case.sn_dec), args.band, args.size, sia, session, cache, limiter))
    with ThreadPoolExecutor(args.workers) as pool:
        futs = {c.case_id: pool.submit(job, c) for c in cases.itertuples()}
        for k, f in futs.items():
            try:
                images[k] = f.result()
            except Exception as e:  # noqa: BLE001 - one bad case must not kill the gallery
                errors[k] = f"{type(e).__name__}: {e}"
                print(f"{k}: {errors[k]}", file=sys.stderr)
    print(f"{len(images)}/{len(cases)} cutouts", file=sys.stderr)
    G.render(cases, images, errors, args)


if __name__ == "__main__":
    main()
