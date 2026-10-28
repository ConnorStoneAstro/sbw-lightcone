"""Truth-to-object match for the OU2024 true hosts (the OU2024 stand-in for
DP0.2's MatchesTruth, which does not exist for OU2024).

For each host (skyCatalog position and r magnitude from ``ou24_sn_table``), every
primary object within ``--match-radius`` arcsec is a candidate. Two matchers:

* ``nearest``: the closest candidate.
* ``mag``: among candidates with |m_r(CModel) - m_r(truth)| < ``--max-dmag``,
  the one with the lowest chi2 = (sep / pos_sigma)^2 + (dmag / mag_sigma)^2.

Flags: n_within (candidates within the radius), multiple (n_within > 1),
no_match (none). Object r magnitudes are AB from cModelFlux (nJy).
Also stores the nearest object's r-band moment shape so truth and pipeline
position angles can be compared (ou24_pa_check). ``truth_match: matches_truth``
selects the ``mag`` matcher, ``nearest_1arcsec`` the ``nearest`` one.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

OBJ_COLS = ["objectId", "coord_ra", "coord_dec", "refExtendedness", "r_cModelFlux", "r_bdReD", "r_bdReB",
            "r_ixx", "r_iyy", "r_ixy", "r_ixxPSF", "r_iyyPSF", "r_ixyPSF", "parentObjectId", "deblend_nChild"]
_H = {}


def unit_vectors(ra_deg, dec_deg):
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _pairs(path):
    t = pq.read_table(path, columns=OBJ_COLS, filters=[("detect_isPrimary", "==", True)]).to_pandas(ignore_metadata=True)
    ra, dec = t.coord_ra.to_numpy(), t.coord_dec.to_numpy()
    h = _H["hosts"]
    pad = _H["radius"] / 3600 + 0.01
    box = np.flatnonzero((h.host_dec_cat >= dec.min() - pad) & (h.host_dec_cat <= dec.max() + pad)
                         & (h.host_ra_cat >= ra.min() - 0.1) & (h.host_ra_cat <= ra.max() + 0.1))
    if len(box) == 0:
        return None
    tree = cKDTree(unit_vectors(ra, dec))
    hu = unit_vectors(h.host_ra_cat.to_numpy()[box], h.host_dec_cat.to_numpy()[box])
    hits = tree.query_ball_point(hu, r=2 * np.sin(np.deg2rad(_H["radius"] / 3600) / 2))
    n = np.fromiter((len(x) for x in hits), int, len(hits))
    if n.sum() == 0:
        return None
    oi = np.concatenate([np.asarray(x, int) for x in hits if len(x)])
    hi = np.repeat(box, n)
    sep = np.rad2deg(2 * np.arcsin(np.linalg.norm(unit_vectors(ra[oi], dec[oi]) - unit_vectors(
        h.host_ra_cat.to_numpy()[hi], h.host_dec_cat.to_numpy()[hi]), axis=1) / 2)) * 3600
    out = t.iloc[oi].reset_index(drop=True)
    out.insert(0, "host_id", h.host_id.to_numpy()[hi])
    out["sep_arcsec"] = sep
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sn-table", required=True, help="ou24_sn_table parquet (host_id, host_ra, host_dec, host_mag_r)")
    ap.add_argument("--objects", required=True, help="objectTable index csv")
    ap.add_argument("--match-radius", type=float, default=1.0, help="arcsec")
    ap.add_argument("--max-dmag", type=float, default=1.0)
    ap.add_argument("--pos-sigma", type=float, default=0.1, help="arcsec, chi2 position scale")
    ap.add_argument("--mag-sigma", type=float, default=0.2, help="mag, chi2 magnitude scale")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    h = pq.read_table(a.sn_table, columns=["host_id", "host_ra", "host_dec", "host_mag_r"]).to_pandas(
        ignore_metadata=True, types_mapper={pa.int64(): pd.Int64Dtype()}.get)     # ids must not pass through float64
    h = h.dropna(subset=["host_id", "host_ra"]).drop_duplicates("host_id")
    h = h.rename(columns={"host_ra": "host_ra_cat", "host_dec": "host_dec_cat"}).reset_index(drop=True)
    h["host_id"] = h.host_id.astype(np.int64)
    _H.update(hosts=h, radius=a.match_radius)
    files = pd.read_csv(a.objects)["path"].tolist()
    with Pool(a.workers) as pool:
        parts = [p for p in pool.map(_pairs, files) if p is not None]
    c = pd.concat(parts, ignore_index=True).drop_duplicates(["host_id", "objectId"])
    c = c.merge(h[["host_id", "host_mag_r"]], on="host_id")
    for k in ("objectId", "parentObjectId"):      # nullable before any left merge: a float round trip corrupts the ids
        c[k] = c[k].astype("Int64")
    with np.errstate(divide="ignore", invalid="ignore"):
        c["obj_mag_r"] = 31.4 - 2.5 * np.log10(c.r_cModelFlux)
    c["dmag"] = c.obj_mag_r - c.host_mag_r
    c["chi2"] = (c.sep_arcsec / a.pos_sigma) ** 2 + (c.dmag / a.mag_sigma) ** 2
    print(f"{len(c)} host-object pairs within {a.match_radius}\" ({time.time() - t0:.0f}s)", file=sys.stderr)

    keep = ["host_id", "objectId", "sep_arcsec", "dmag", "obj_mag_r", "refExtendedness", "r_bdReD", "r_bdReB",
            "r_ixx", "r_iyy", "r_ixy", "r_ixxPSF", "r_iyyPSF", "r_ixyPSF", "parentObjectId", "deblend_nChild"]
    near = c.sort_values("sep_arcsec").drop_duplicates("host_id")[keep]
    near = near.rename(columns={k: f"nearest_{k}" for k in keep if k != "host_id"})
    mm = c[np.abs(c.dmag) < a.max_dmag].sort_values("chi2").drop_duplicates("host_id")
    mm = mm[["host_id", "objectId", "sep_arcsec", "dmag", "chi2"]].rename(
        columns={"objectId": "mag_objectId", "sep_arcsec": "mag_sep_arcsec", "dmag": "mag_dmag", "chi2": "mag_chi2"})
    nw = c.groupby("host_id").size().rename("n_within").reset_index()
    out = h.merge(nw, on="host_id", how="left").merge(near, on="host_id", how="left").merge(mm, on="host_id", how="left")
    out["n_within"] = out.n_within.fillna(0).astype(int)
    out["no_match"] = out.n_within == 0
    out["multiple"] = out.n_within > 1
    assert all(out[k].dtype == "Int64" for k in ("nearest_objectId", "mag_objectId", "nearest_parentObjectId"))
    out["matchers_agree"] = out.nearest_objectId.fillna(-1).to_numpy() == out.mag_objectId.fillna(-1).to_numpy()
    out["host_id"] = out.host_id.astype("Int64")
    tab = pa.Table.from_pandas(out, preserve_index=False).replace_schema_metadata(
        {b"match_radius": str(a.match_radius).encode(), b"max_dmag": str(a.max_dmag).encode(),
         b"pos_sigma": str(a.pos_sigma).encode(), b"mag_sigma": str(a.mag_sigma).encode()})
    pq.write_table(tab, a.output, compression="zstd")
    n = len(out)
    hasn, hasm = out.nearest_objectId.notna(), out.mag_objectId.notna()
    print(f"hosts {n}: nearest match {hasn.mean():.4f}, mag match {hasm.mean():.4f}, multiple {out.multiple.mean():.4f}; "
          f"both matched but disagree {(hasn & hasm & ~out.matchers_agree).sum()} "
          f"({(hasn & hasm & ~out.matchers_agree).sum() / max((hasn & hasm).sum(), 1):.4f}); "
          f"nearest only {(hasn & ~hasm).sum()} ({time.time() - t0:.0f}s)", file=sys.stderr)


if __name__ == "__main__":
    main()
