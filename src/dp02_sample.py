"""Build the DP0.2 SN test sample: SN truth, true host truth, and the host's
DP0.2 objectId, with the sample cuts applied."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

FAR_ARCSEC = 3.0


def ang_sep(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.deg2rad, (ra1, dec1, ra2, dec2))
    h = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
    return np.rad2deg(2 * np.arcsin(np.sqrt(h))) * 3600


def _sne(f):
    t = pq.read_table(f, columns=["id", "id_string", "host_galaxy", "truth_type", "ra", "dec", "redshift", "tract"],
                      filters=[("truth_type", "==", 3)])
    return t.drop_columns(["truth_type"]).to_pandas(ignore_metadata=True)


def _hosts(args):
    f, ids = args
    t = pq.read_table(f, columns=["id", "truth_type", "ra", "dec", "redshift", "mag_r", "flux_r"],
                      filters=[("truth_type", "==", 1)])
    t = t.filter(pc.is_in(t["id"], value_set=pa.array(ids, pa.int64())))
    return t.drop_columns(["truth_type"]).to_pandas(ignore_metadata=True)


def _match(args):
    f, idstr = args
    t = pq.read_table(f, columns=["id", "match_objectId", "match_candidate", "match_count", "match_chisq"])
    t = t.filter(pc.is_in(t["id"], value_set=pa.array(idstr, pa.string())))
    return t.to_pandas(ignore_metadata=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--match", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--truth-match", required=True, choices=["matches_truth", "nearest_1arcsec"])
    ap.add_argument("--cuts", required=True, choices=["clean", "all"])
    ap.add_argument("--workers", type=int, default=int(os.environ.get("DP02_WORKERS", "4")))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    tfiles = sorted(glob.glob(os.path.join(a.truth, "truth_tract*.parquet")))
    with Pool(a.workers) as pool:
        sn = pd.concat(pool.map(_sne, tfiles), ignore_index=True).drop_duplicates("id")
        print(f"{len(sn)} SNe, {(sn.host_galaxy > 0).sum()} hosted, {time.time() - t0:.0f}s", file=sys.stderr)
        hid = np.unique(sn.loc[sn.host_galaxy > 0, "host_galaxy"].to_numpy())
        hosts = pd.concat(pool.map(_hosts, [(f, hid) for f in tfiles]), ignore_index=True).drop_duplicates("id")
        print(f"{len(hosts)} host rows for {len(hid)} host ids, {time.time() - t0:.0f}s", file=sys.stderr)
        mfiles = sorted(glob.glob(os.path.join(a.match, "*.parq")))
        idstr = hid.astype(str)
        m = pd.concat(pool.map(_match, [(f, idstr) for f in mfiles]), ignore_index=True)
    m = m[(m.match_objectId > 0)].sort_values("match_chisq").drop_duplicates("id")
    m["host_galaxy"] = m["id"].astype(np.int64)
    print(f"{len(m)} hosts matched to DP0.2 objects, {time.time() - t0:.0f}s", file=sys.stderr)

    hosts = hosts.rename(columns={"id": "host_galaxy", "ra": "host_ra", "dec": "host_dec",
                                  "redshift": "host_redshift", "mag_r": "host_mag_r", "flux_r": "host_flux_r"})
    df = sn.rename(columns={"id": "sn_id", "id_string": "sn_id_string", "ra": "sn_ra", "dec": "sn_dec",
                            "redshift": "sn_redshift"})
    df = df.merge(hosts, on="host_galaxy", how="left")
    df["hostless"] = df.host_galaxy <= 0
    df["host_missing"] = (~df.hostless) & (df.host_ra.isna() | ~(df.host_flux_r > 0))
    df["sep_truth_arcsec"] = ang_sep(df.sn_ra, df.sn_dec, df.host_ra, df.host_dec)
    df["far"] = df.sep_truth_arcsec > FAR_ARCSEC

    # both truth matchers are always computed; host_objectId is the one the decision selects
    mm = m[["host_galaxy", "match_objectId", "match_chisq"]].copy()
    mm["match_objectId"] = mm["match_objectId"].astype("Int64")   # nullable: a float merge would round the ids
    df = df.merge(mm, on="host_galaxy", how="left")
    df["host_objectId_mt"] = df.match_objectId.fillna(-1).astype(np.int64)
    df = df.drop(columns=["match_objectId"])
    c = pq.read_table(a.candidates, columns=["sn_id", "objectId", "coord_ra", "coord_dec"]).to_pandas(ignore_metadata=True)
    c = c.merge(df[["sn_id", "host_ra", "host_dec"]], on="sn_id")
    c["dh"] = ang_sep(c.coord_ra, c.coord_dec, c.host_ra, c.host_dec)
    c = c[c.dh <= 1.0].sort_values("dh").drop_duplicates("sn_id")
    c["objectId"] = c["objectId"].astype("Int64")
    df = df.merge(c[["sn_id", "objectId", "dh"]].rename(columns={"objectId": "host_objectId_nn", "dh": "match_sep_arcsec"}),
                  on="sn_id", how="left")
    df["host_objectId_nn"] = df.host_objectId_nn.fillna(-1).astype(np.int64)
    df["host_objectId"] = df["host_objectId_mt" if a.truth_match == "matches_truth" else "host_objectId_nn"]
    df["host_unmatched"] = (~df.hostless) & (df.host_objectId < 0)

    # how the two matchers compare on hosted SNe with a host row, before the matcher-dependent cut
    e = df[~df.hostless & ~df.host_missing & (~df.far if a.cuts != "all" else True)]
    mt, nn = e.host_objectId_mt.to_numpy() >= 0, e.host_objectId_nn.to_numpy() >= 0
    precut = {"n": len(e), "both": int((mt & nn).sum()), "mt_only": int((mt & ~nn).sum()),
              "nn_only": int((~mt & nn).sum()), "neither": int((~mt & ~nn).sum()),
              "both_differ": int((mt & nn & (e.host_objectId_mt != e.host_objectId_nn).to_numpy()).sum())}
    print(f"truth matchers before cuts: {precut}", file=sys.stderr)

    n0 = len(df)
    if a.cuts == "clean":
        df = df[~df.hostless & ~df.host_missing & ~df.host_unmatched & ~df.far]
    df = df.sort_values("sn_id").reset_index(drop=True)
    print(f"sample: {len(df)} of {n0} SNe after cuts={a.cuts} ({time.time() - t0:.0f}s)", file=sys.stderr)
    tab = pa.Table.from_pandas(df, preserve_index=False)
    tab = tab.replace_schema_metadata({**tab.schema.metadata, b"truth_match_precut": json.dumps(precut).encode()})
    pq.write_table(tab, a.output, compression="zstd")


if __name__ == "__main__":
    main()
