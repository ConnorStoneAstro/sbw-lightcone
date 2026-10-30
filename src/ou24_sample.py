"""Build the OU2024 transient test sample: transient truth, true host truth and the host's objectId
under truth_match, with the dp02_sample_cuts options. Output columns follow dp02_sample, so the DP0.2
association, sweep and figure scripts run on it unchanged.

truth_match (OU2024 has no MatchesTruth table; both come from ou24_truth_match):
  matches_truth   -> the magnitude-aware match (|dm_r| < 1 within 1 arcsec, lowest chi2)
  nearest_1arcsec -> the nearest primary object within 1 arcsec
Both are always carried (host_objectId_mt, host_objectId_nn).

Offsets: offset_re_bulge / offset_re_disk are the SN's elliptical radius around the host truth position
in units of the skyCatalog bulge / disk half-light semi-major axis (semi-major = hlr / sqrt(q)), with
the major axis at host_beta (from +x West toward +y North; SNANA placed transients on this ellipse, see
ou24_pa_check). offset_re uses the dominant component in i (B/T_i > 0.5: bulge, else disk), as
dp02_sample does with cosmoDC2 B/T_i. FIXMAG transients peak after the survey and are always dropped.

As dp02_sample, the sample is rung 2 of the ladder: hosted transients with a host row, detectable
(transient_peak_mag on the SNANA PEAKMAG_r; every OU2024 transient peaks inside the survey) and within the
offset cut; host-unmatched transients stay with host_objectId = -1; metadata key sample_ladder has the counts.
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from dp02_lib import BULGE_DOMINATED_BT, PEAK_MAG
from dp02_sample import FAR_ARCSEC, N_RE, ang_sep
from sbw_host.core import sky_to_plane

INT64 = {pa.int64(): pd.Int64Dtype()}.get     # nullable int64 columns must not come back as float64 (19-digit ids)
COLS = ["sn_id", "sn_ra", "sn_dec", "sn_redshift", "class", "tract", "host_id", "hostless", "host_missing",
        "GALSNSEP", "GALSNDDLR", "a_rot", "peak_mjd", "PEAKMAG_r", "host_ra", "host_dec", "host_redshift",
        "host_mag_r", "host_mag_i", "host_hlr_bulge", "host_hlr_disk", "host_q_bulge", "host_q_disk",
        "host_beta_bulge", "host_beta_disk", "host_bt_r", "host_bt_i", "host_logmass_cat"]


def host_offsets(df):
    """Elliptical SN offset in units of the bulge, disk and dominant-component half-light semi-major axis."""
    x, y = sky_to_plane(df["sn_ra"].to_numpy(float), df["sn_dec"].to_numpy(float),
                        df["host_ra"].to_numpy(float), df["host_dec"].to_numpy(float))   # arcsec, +x West, +y North
    out = {}
    for comp in ("bulge", "disk"):
        th = np.deg2rad(df[f"host_beta_{comp}"].to_numpy(float))
        q = df[f"host_q_{comp}"].to_numpy(float)
        hlr = df[f"host_hlr_{comp}"].to_numpy(float)
        u, v = x * np.cos(th) + y * np.sin(th), -x * np.sin(th) + y * np.cos(th)
        with np.errstate(invalid="ignore", divide="ignore"):
            a = hlr / np.sqrt(q)
            out[comp] = np.where(a > 0, np.hypot(u, v / q) / a, np.nan)
    bd = df["host_bt_i"].to_numpy(float) > BULGE_DOMINATED_BT
    out["own"] = np.where(bd, out["bulge"], out["disk"])
    out["bulge_dominated"] = bd
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sn-table", required=True)
    ap.add_argument("--match", required=True, help="ou24_truth_match")
    ap.add_argument("--truth-match", required=True, choices=["matches_truth", "nearest_1arcsec"])
    ap.add_argument("--cuts", required=True, choices=["clean3", "clean10", "physical", "all"])
    ap.add_argument("--peak-mag", required=True, choices=sorted(PEAK_MAG))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = pq.read_table(a.sn_table, columns=COLS).to_pandas(ignore_metadata=True, types_mapper=INT64)
    m = pq.read_table(a.match, columns=["host_id", "nearest_objectId", "nearest_sep_arcsec", "mag_objectId",
                                        "n_within", "multiple"]).to_pandas(ignore_metadata=True, types_mapper=INT64)
    assert m.nearest_objectId.dtype == "Int64" and m.mag_objectId.dtype == "Int64"
    n_truth = len(df)
    df = df[df["class"] != "FIXMAG"]
    df = df.merge(m, on="host_id", how="left")
    df["host_objectId_mt"] = df.mag_objectId.fillna(-1).astype(np.int64)
    df["host_objectId_nn"] = df.nearest_objectId.fillna(-1).astype(np.int64)
    df = df.drop(columns=["mag_objectId", "nearest_objectId"]).rename(columns={"nearest_sep_arcsec": "match_sep_arcsec"})
    df["host_objectId"] = df["host_objectId_mt" if a.truth_match == "matches_truth" else "host_objectId_nn"]
    df["host_galaxy"] = df.host_id
    df["host_unmatched"] = ~df.hostless & (df.host_objectId < 0)
    df["sep_truth_arcsec"] = ang_sep(df.sn_ra, df.sn_dec, df.host_ra, df.host_dec)
    df["peak_mag_r"] = df.PEAKMAG_r.to_numpy(float)
    lim = PEAK_MAG[a.peak_mag]
    df["faint"] = False if lim is None else ~(df.peak_mag_r < lim).to_numpy()
    off = host_offsets(df)
    df["offset_re_bulge"], df["offset_re_disk"], df["offset_re"] = off["bulge"], off["disk"], off["own"]
    df["bulge_dominated"] = off["bulge_dominated"]
    if a.cuts in FAR_ARCSEC:
        df["far"] = df.sep_truth_arcsec > FAR_ARCSEC[a.cuts]
    elif a.cuts == "physical":
        df["far"] = df.offset_re.to_numpy() > N_RE
        print(f"physical cut N = {N_RE:g}: {int(df.far.sum())} transients beyond, "
              f"{int((~df.hostless & ~df.host_missing & df.offset_re.isna()).sum())} hosted without shapes kept",
              file=sys.stderr)
    else:
        df["far"] = False

    hosted = ~df.hostless
    row = hosted & ~df.host_missing
    det = row & ~df.faint
    samp = det & ~df.far
    ladder = {"transients": int(len(df)), "transients_with_fixmag": int(n_truth), "hosted": int(hosted.sum()),
              "host_row": int(row.sum()), "detectable": int(det.sum()), "sample": int(samp.sum()),
              "faint_dropped": int((row & df.faint).sum()), "far_dropped": int((det & df.far).sum()),
              "host_unmatched_in_sample": int((samp & df.host_unmatched).sum()),
              "peak_mag_cut": a.peak_mag, "peak_mag_limit": lim, "cuts": a.cuts,
              "n_re": N_RE if a.cuts == "physical" else None}
    print(f"ladder: {ladder}", file=sys.stderr)
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
    df["sn_id"] = df.sn_id.astype(np.int64)
    print(f"sample: {len(df)} of {n0} transients (FIXMAG dropped; {n_truth} in the table) after cuts={a.cuts}, "
          f"peak {a.peak_mag}", file=sys.stderr)
    tab = pa.Table.from_pandas(df, preserve_index=False)
    tab = tab.replace_schema_metadata({**tab.schema.metadata, b"truth_match_precut": json.dumps(precut).encode(),
                                       b"sample_ladder": json.dumps(ladder).encode()})
    pq.write_table(tab, a.output, compression="zstd")


if __name__ == "__main__":
    main()
