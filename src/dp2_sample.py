"""The DP2 transient sample chosen by decision dp2_sample, as a table keyed like the DP0.2 sample.

Reads the EDP2 transient working-group file named by the option from the index of those files
(data/dp2_sn_candidates_index.csv), after checking that its size and mtime still match the index, so a
re-release by the collaborators fails loudly instead of silently changing the sample. Writes sn_id (the
DiaObject id), sn_ra, sn_dec, host_objectId = -1 (no truth) and the light-curve fit columns the file has
(SALT3 z, x1, c, mB, FITPROB, ... for rc_salt).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

FILES = {"rc_salt": "RC_sample.parquet", "md_final": "final_sample.parquet", "parsnip_ddf": "parsnip_ddf.parquet",
         "rc_all": "RC_lcfit_input_122k.parquet"}
# (id, ra, dec) column names per file
POS = {"rc_salt": ("diaObjectId", "ra", "dec"),
       "md_final": ("diaObjectId_dia_object_lc", "ra_dia_object_lc", "dec_dia_object_lc"),
       "parsnip_ddf": ("diaObjectId", "ra", "dec"),
       "rc_all": ("diaObjectId", "ra", "dec")}
KEEP = ["z", "z_err", "PKMJD", "x1", "x1ERR", "c", "cERR", "mB", "mBERR", "FITPROB", "MWEBV", "SNRMAX1", "SNRSUM",
        "NOBS_DETECT", "N_NIGHTS_SNR5", "N_BANDS_SNR5", "in_MD", "MD_set", "MD_parsnip_p_SNIa",
        "predicted_redshift", "count"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", required=True, help="csv (path, size, mtime) of the EDP2 sample files")
    ap.add_argument("--sample", required=True, choices=sorted(FILES))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    idx = pd.read_csv(a.index)
    row = idx[idx.path.str.endswith("/" + FILES[a.sample])]
    if len(row) != 1:
        raise SystemExit(f"{FILES[a.sample]} is not in {a.index}")
    path, size, mtime = row.iloc[0][["path", "size", "mtime"]]
    st = os.stat(path)
    if st.st_size != size or int(st.st_mtime) != mtime:
        raise SystemExit(f"{path} changed since the index was written (size {st.st_size} vs {size}, "
                         f"mtime {int(st.st_mtime)} vs {mtime}); rewrite the index deliberately")
    sid, ra, dec = POS[a.sample]
    names = pq.read_schema(path).names
    df = pq.read_table(path, columns=[sid, ra, dec] + [c for c in KEEP if c in names]).to_pandas(ignore_metadata=True)
    df = df.rename(columns={sid: "sn_id", ra: "sn_ra", dec: "sn_dec"})
    df["sn_id"] = df["sn_id"].astype(np.int64)
    if df.sn_id.duplicated().any():
        raise SystemExit(f"duplicate DiaObject ids in {path}")
    df["host_objectId"] = np.int64(-1)
    df = df.sort_values("sn_id").reset_index(drop=True)
    df.attrs = {"dp2_sample": a.sample, "source": path}
    df.to_parquet(a.output, index=False)
    print(f"{a.sample}: {len(df)} transients from {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
