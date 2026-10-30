"""Build the DP0.2 SN test sample: SN truth, true host truth, and the host's
DP0.2 objectId, with the sample cuts applied.

The sample (rung 2 of the ladder) is every hosted SN with a host row that is detectable (transient_peak_mag:
peak r from the DC2 SN truth summary, with the peak epoch inside the DP0.2 visit window) and within the
dp02_sample_cuts offset cut. SNe whose true host has no DP0.2 object (host_unmatched) stay in the sample
with host_objectId = -1, so the association records what the pipeline does with them. The parquet
metadata key sample_ladder carries the counts at every step, including the SNe dropped as too faint."""
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

from dp02_lib import PEAK_MAG, host_offsets

FAR_ARCSEC = {"clean3": 3.0, "clean10": 10.0}   # geometric cuts: drop SNe this far from the truth host
# physical cut: drop SNe beyond N_RE half-light radii of the host along its truth ellipse (bulge R_e for
# bulge-dominated hosts, disk R_e otherwise; see dp02_offset_diagnostic). Hosts without cosmoDC2 shapes are kept.
N_RE = 3.0   # where both-wrong-same-object turns up in dp02_offset_diagnostic (sep10): ~3% at 2 R_e, ~12% at 3.6 R_e
SHAPE_COLS = ["galaxy_id", "size_bulge_true", "size_minor_bulge_true", "size_disk_true", "size_minor_disk_true",
              "position_angle_true", "bulge_to_total_ratio_i"]


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
    ap.add_argument("--shapes", required=True, help="dp02_host_truth_shapes (used by cuts=physical)")
    ap.add_argument("--peak", required=True, help="DC2 SN truth summary parquet (t0, max_flux_r per SN)")
    ap.add_argument("--visits", required=True, help="DP0.2 visitTable parquet (expMidptMJD): the imaging window")
    ap.add_argument("--peak-mag", required=True, choices=sorted(PEAK_MAG))
    ap.add_argument("--truth-match", required=True, choices=["matches_truth", "nearest_1arcsec"])
    ap.add_argument("--cuts", required=True, choices=["clean3", "clean10", "physical", "all"])
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
    # detectability: true peak r (10-year maximum flux) and the peak epoch inside the DP0.2 visit window
    pk = pq.read_table(a.peak, columns=["id", "t0", "max_flux_r"]).to_pandas(ignore_metadata=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        pk["peak_mag_r"] = np.where(pk.max_flux_r > 0, 31.4 - 2.5 * np.log10(pk.max_flux_r), np.nan)
    pk = pk.rename(columns={"id": "sn_id", "t0": "peak_mjd"})[["sn_id", "peak_mjd", "peak_mag_r"]]
    df = df.merge(pk, on="sn_id", how="left")
    mjd = pq.read_table(a.visits, columns=["expMidptMJD"]).column(0).to_numpy()
    window = (float(mjd.min()), float(mjd.max()))
    df["in_window"] = (df.peak_mjd >= window[0]) & (df.peak_mjd <= window[1])
    lim = PEAK_MAG[a.peak_mag]
    df["faint"] = False if lim is None else ~(df.in_window & (df.peak_mag_r < lim)).to_numpy()
    if a.cuts in FAR_ARCSEC:
        df["far"] = df.sep_truth_arcsec > FAR_ARCSEC[a.cuts]
    elif a.cuts == "physical":
        sh = pq.read_table(a.shapes, columns=SHAPE_COLS).to_pandas(ignore_metadata=True)
        d = df[["host_galaxy", "sn_ra", "sn_dec", "host_ra", "host_dec"]].merge(
            sh.rename(columns={"galaxy_id": "host_galaxy"}), on="host_galaxy", how="left")
        df["offset_re"] = host_offsets(d)["own"]
        df["far"] = df.offset_re.to_numpy() > N_RE
        print(f"physical cut N = {N_RE:g}: {int(df.far.sum())} SNe beyond, "
              f"{int((~df.hostless & ~df.host_missing & df.offset_re.isna()).sum())} hosted without shapes kept", file=sys.stderr)
    else:
        df["far"] = False

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

    hosted = ~df.hostless
    row = hosted & ~df.host_missing
    det = row & ~df.faint
    samp = det & ~df.far
    ladder = {"transients": int(len(df)), "hosted": int(hosted.sum()), "host_row": int(row.sum()),
              "detectable": int(det.sum()), "sample": int(samp.sum()),
              "faint_dropped": int((row & df.faint).sum()), "far_dropped": int((det & df.far).sum()),
              "host_unmatched_in_sample": int((samp & df.host_unmatched).sum()),
              "peak_mag_cut": a.peak_mag, "peak_mag_limit": lim, "visit_window_mjd": list(window),
              "cuts": a.cuts, "n_re": N_RE if a.cuts == "physical" else None}
    print(f"ladder: {ladder}", file=sys.stderr)
    # how the two matchers compare on the sample, before the matcher-dependent split
    e = df[samp]
    mt, nn = e.host_objectId_mt.to_numpy() >= 0, e.host_objectId_nn.to_numpy() >= 0
    precut = {"n": len(e), "both": int((mt & nn).sum()), "mt_only": int((mt & ~nn).sum()),
              "nn_only": int((~mt & nn).sum()), "neither": int((~mt & ~nn).sum()),
              "both_differ": int((mt & nn & (e.host_objectId_mt != e.host_objectId_nn).to_numpy()).sum())}
    print(f"truth matchers before cuts: {precut}", file=sys.stderr)

    n0 = len(df)
    if a.cuts != "all":
        df = df[samp]
    df = df.sort_values("sn_id").reset_index(drop=True)
    print(f"sample: {len(df)} of {n0} SNe after cuts={a.cuts}, peak {a.peak_mag} ({time.time() - t0:.0f}s)", file=sys.stderr)
    tab = pa.Table.from_pandas(df, preserve_index=False)
    tab = tab.replace_schema_metadata({**tab.schema.metadata, b"truth_match_precut": json.dumps(precut).encode(),
                                       b"sample_ladder": json.dumps(ladder).encode()})
    pq.write_table(tab, a.output, compression="zstd")


if __name__ == "__main__":
    main()
