"""One-off fetch of every DES DR2 object within 2 arcmin of every DES-SN5YR Dovekie SN (IDSURVEY = 10) from
NOIRLab Data Lab TAP (des_dr2.main), written once to data/des_dr2_candidates.parquet. Not a recipe: TAP is
not reproducible, so the cache is the declared input des_dr2_objects and its fetch is recorded in
data/README.md and in the parquet metadata (query template, date, service, SN count).

One row per unique coadd_object_id (objects shared by overlapping cones are kept once); src/des_candidates.py
recomputes cone membership. 2 arcmin covers the largest search radius used for the SNe (60") and random
positions drawn 45-90" from an SN with a 30" cone.

    python src/des_fetch_dr2.py --metadata <DES-Dovekie_Metadata.csv> --head <DES-SN5YR_DES_HEAD.FITS.gz> \
        --output data/des_dr2_candidates.parquet
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyvo

TAP = "https://datalab.noirlab.edu/tap"
RADIUS_DEG = 2.0 / 60
BANDS = "grizy"
GENERAL = ["coadd_object_id", "ra", "dec", "alphawin_j2000", "deltawin_j2000", "a_image", "b_image",
           "theta_j2000", "kron_radius", "extended_class_coadd", "extended_class_wavg", "ebv_sfd98"]
PER_BAND = ["flux_auto", "fluxerr_auto", "mag_auto", "flux_radius", "kron_radius", "awin_image", "bwin_image",
            "thetawin_image", "class_star", "spread_model", "spreaderr_model", "imaflags_iso", "flags", "niter_model"]
COLUMNS = GENERAL + [f"{c}_{b}" for b in BANDS for c in PER_BAND]
QUERY = "SELECT {cols} FROM des_dr2.main WHERE q3c_radial_query(ra, dec, {ra:.7f}, {dec:.7f}, {r:.7f}) = 't'"


def sn_positions(metadata, head):
    """Dovekie DES SNe (IDSURVEY = 10) with the HEAD-file SN position."""
    from astropy.table import Table
    m = pd.read_csv(metadata, sep=r"\s+", comment="#")
    m = m[(m["VARNAMES:"] == "SN:") & (m.IDSURVEY.astype(int) == 10)]
    h = Table.read(head)
    h = h[["SNID", "RA", "DEC"]].to_pandas()
    h["SNID"] = h.SNID.apply(lambda s: s.decode() if isinstance(s, bytes) else s).str.strip().astype(np.int64)
    out = m[["CID"]].astype(np.int64).merge(h, left_on="CID", right_on="SNID", how="left")
    if out.RA.isna().any():
        raise SystemExit(f"{out.RA.isna().sum()} Dovekie SNe missing from the HEAD file")
    return out


def cone(svc, ra, dec, tries=6):
    q = QUERY.format(cols=", ".join(COLUMNS), ra=ra, dec=dec, r=RADIUS_DEG)
    for k in range(tries):
        try:
            return svc.search(q, maxrec=100000).to_table().to_pandas()
        except Exception as e:   # transient TAP errors: back off and retry
            if k == tries - 1:
                raise
            print(f"retry {k + 1} ({ra:.4f}, {dec:.4f}): {e}", file=sys.stderr)
            time.sleep(5 * (k + 1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--head", required=True)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    sn = sn_positions(a.metadata, a.head)
    svc = pyvo.dal.TAPService(TAP)
    t0 = time.time()
    parts, nper = [], []
    with ThreadPoolExecutor(a.threads) as ex:
        for i, df in enumerate(ex.map(lambda r: cone(svc, r.RA, r.DEC), sn.itertuples())):
            parts.append(df); nper.append(len(df))
            if (i + 1) % 100 == 0:
                print(f"{i + 1}/{len(sn)} cones, {time.time() - t0:.0f}s", file=sys.stderr)
    nper = np.array(nper)
    if (nper >= 100000).any():
        raise SystemExit("a cone hit maxrec")
    obj = pd.concat(parts, ignore_index=True).drop_duplicates("coadd_object_id").reset_index(drop=True)
    obj["coadd_object_id"] = obj.coadd_object_id.astype(np.int64)
    meta = {"service": TAP, "table": "des_dr2.main", "query_template": QUERY.format(cols="<columns>", ra=0, dec=0, r=RADIUS_DEG),
            "radius_arcsec": str(RADIUS_DEG * 3600), "n_sn": str(len(sn)), "n_cone_rows": str(int(nper.sum())),
            "n_objects": str(len(obj)), "min_per_cone": str(int(nper.min())), "median_per_cone": str(float(np.median(nper))),
            "fetched_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "positions": "DES-SN5YR_DES_HEAD.FITS.gz RA, DEC of the Dovekie IDSURVEY = 10 SNe"}
    tab = pa.Table.from_pandas(obj, preserve_index=False)
    tab = tab.replace_schema_metadata({k.encode(): v.encode() for k, v in meta.items()})
    pq.write_table(tab, a.output, compression="zstd")
    print(meta, f"{time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
