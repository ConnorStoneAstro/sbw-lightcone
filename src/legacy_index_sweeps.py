"""Write the index (path, size, mtime, sha256) of the Legacy Surveys DR10 south sweep tiles that hold any
position of the given parquet files (sn_ra, sn_dec) padded by --pad-deg, which astra.yaml declares as the
input legacy_dr10_sweeps in place of the 1437-tile directory. Each sha256 is checked against the
survey's own legacysurvey_dr10_south_sweep_10.1.sha256sum in the same directory.

    python src/legacy_index_sweeps.py --sweep-dir <dir> --pad-deg 0.3 --output data/legacy_dr10_sweeps_index.csv \
        results/dp2_ddf/dp2_sample.parquet [...]

Sweep tiles are 5 x 5 deg boxes named sweep-<ramin><p|m><decmin>-<ramax><p|m><decmax>.fits
(src/dp2_candidates_legacy.py uses the same naming to pick the tiles of a position).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys

import numpy as np
import pandas as pd


def tile_name(ra0, dec0):
    """Sweep file name of the 5 x 5 deg tile with lower corner (ra0, dec0)."""
    s = lambda d: f"{'p' if d >= 0 else 'm'}{abs(int(d)):03d}"
    return f"sweep-{int(ra0):03d}{s(dec0)}-{int(ra0) + 5:03d}{s(dec0 + 5)}.fits"


def tiles_near(ra, dec, pad_deg):
    """Names of every tile within pad_deg of (ra, dec) (RA padding scaled by 1/cos dec)."""
    padra = pad_deg / max(np.cos(np.deg2rad(min(abs(dec) + pad_deg, 89.9))), 1e-3)
    out = set()
    for r in (ra - padra, ra, ra + padra):
        for d in (dec - pad_deg, dec, dec + pad_deg):
            out.add(tile_name(5 * np.floor((r % 360) / 5), 5 * np.floor(np.clip(d, -90, 89.999) / 5)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep-dir", required=True)
    ap.add_argument("--pad-deg", type=float, default=0.3)
    ap.add_argument("--output", required=True)
    ap.add_argument("positions", nargs="+")
    a = ap.parse_args(argv)
    pos = pd.concat([pd.read_parquet(p, columns=["sn_ra", "sn_dec"]) for p in a.positions], ignore_index=True)
    names = sorted(set().union(*[tiles_near(r, d, a.pad_deg) for r, d in zip(pos.sn_ra, pos.sn_dec)]))
    ref = {}
    for line in open(os.path.join(a.sweep_dir, "legacysurvey_dr10_south_sweep_10.1.sha256sum")):
        h, f = line.split()
        ref[os.path.basename(f.lstrip("*"))] = h
    rows, absent = [], []
    for n in names:
        p = os.path.join(a.sweep_dir, n)
        if not os.path.exists(p):
            absent.append(n)
            continue
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for blk in iter(lambda: f.read(1 << 24), b""):
                h.update(blk)
        if ref.get(n) != h.hexdigest():
            raise SystemExit(f"{n}: sha256 differs from the survey's sha256sum")
        st = os.stat(p)
        rows.append([p, st.st_size, int(st.st_mtime), h.hexdigest()])
    with open(a.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["path", "size", "mtime", "sha256"])
        w.writerows(rows)
    print(f"{len(pos)} positions: {len(rows)} tiles, {sum(r[1] for r in rows) / 1e9:.1f} GB; "
          f"no sweep file (outside DR10 south) for {absent}", file=sys.stderr)


if __name__ == "__main__":
    main()
