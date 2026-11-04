"""SBW and DLR on DP2 against Y. J. Kang's DELIGHT host predictions (no truth: agreement, never recovery).

DELIGHT (Forster+22) ran on Legacy Surveys / PS1 / DES DR2 cutouts of each transient, not on Rubin pixels, and
always returns a position: the CNN host centre (ra_delight, dec_delight, 8-fold averaged) and the SExtractor
source nearest in elliptical distance (ra_sex, dec_sex, "DELIGHT-SEx"). Both positions are matched to the DP2
candidates of the transient (every primary object within the search radius, before any object-table cut)
under decision reference_match: pos_1arcsec / pos_3arcsec take the nearest candidate within 1" / 3",
frac_semimajor the candidate whose half-light ellipse (sbw_profile, headline band) holds the position at the
smallest elliptical radius <= 0.5 semi-major axes. A position with no candidate is DELIGHT's implicit "no
host". SBW is hostless when flb_best < sbw_hostless_cut, DLR when no candidate has d_DLR <= dlr_max.

The association is recomputed with the DELIGHT-matched object as ``host_objectId``, which gives P_SBW and d_DLR
of that object and whether the object-table cuts removed it (``true_host_status``); the SBW and DLR picks are
checked against dp2_associations. Modes:

* ``agreement``: csv, one row per (image_survey subset, method pair): same host / different / one-hostless /
  both-hostless counts and fractions; subsets ``all`` and each image_survey.
* ``metrics``: json, the headline fractions with Wilson 68% intervals, per image_survey, hostless rates per
  method, the match cross-check against Kang's dp2_objectId, the P_SBW(DELIGHT host) distribution, and the
  SBW / DLR agreement with DELIGHT under every reference_match rule.
* ``cases``: csv for src/dp02_gallery.py, 12 seeded transients per disagreement class with SBW / DLR hosts
  and ellipses and the raw DELIGHT and DELIGHT-SEx positions (markers ``delight``, ``sex``).

Only transients present both in dp2_associations and in Kang's table are used; in a universe whose dp2_sample
is not parsnip_ddf the overlap may be small or empty, and the outputs are then empty but valid.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import dp02_lib as L
import dp2_lib
from dp02_lib import HOSTLESS_CUT, P_TRUE, associate_all, sbw_hostless, wilson
from sbw_host.core import sky_to_plane

MATCH = {"pos_1arcsec": 1.0, "pos_3arcsec": 3.0, "frac_semimajor": 0.5}
METHODS = ("sbw", "dlr", "delight", "sex")
LABEL = {"sbw": "SBW", "dlr": "DLR", "delight": "DELIGHT", "sex": "DELIGHT-SEx"}
PAIRS = [("sbw", "dlr"), ("sbw", "delight"), ("sbw", "sex"), ("dlr", "delight"), ("dlr", "sex"), ("delight", "sex")]
CELLS = ("same", "different", "a_only_hostless", "b_only_hostless", "both_hostless")
N_PER_CLASS = 12
SEED = 20261029
RESULTS, COMPANION = "delight_parsnip_dp2_8k_results.parquet", "delight_parsnip_dp2_8k_legacy_dr10_photoz_specz.parquet"


def read_kang(index, ls_id_shift=20):
    """Kang's two tables, after checking size, mtime and sha256 against the index; one row per transient."""
    idx = pd.read_csv(index)
    paths = {}
    for name in (RESULTS, COMPANION):
        row = idx[idx.path.str.endswith("/" + name)]
        if len(row) != 1:
            raise SystemExit(f"{name} is not in {index}")
        path, size, mtime, sha = row.iloc[0][["path", "size", "mtime", "sha256"]]
        st = os.stat(path)
        if st.st_size != size or int(st.st_mtime) != mtime:
            raise SystemExit(f"{path} changed since the index was written; rewrite the index deliberately")
        if hashlib.sha256(open(path, "rb").read()).hexdigest() != sha:
            raise SystemExit(f"{path} content differs from the index sha256")
        paths[name] = path
    res = pd.read_parquet(paths[RESULTS], columns=["oid", "image_survey", "ra_delight", "dec_delight", "ra_sex",
                                                   "dec_sex", "std_delight", "hostsize"])
    com = pd.read_parquet(paths[COMPANION], columns=["diaobjectId", "dp2_objectId", "match_sep_arcsec", "legacy_ls_id"])
    com = com.rename(columns={"diaobjectId": "oid", "dp2_objectId": "kang_dp2_objectId",
                              "match_sep_arcsec": "kang_match_sep"})
    com["kang_dp2_objectId"] = com.kang_dp2_objectId.fillna(-1).astype(np.int64)
    com["kang_match_sep"] = com.kang_match_sep.astype(float)
    # Kang's legacy_ls_id (a string) = release << 40 | brickid << ls_id_shift | objid; with shift 20 (bits 16-19
    # are always zero) 99.7% decode to a dp2_ddf_legacy candidate at the identical position, with 16 none do.
    # -> the objectId of src/dp2_candidates_legacy.py
    ls = pd.to_numeric(com.pop("legacy_ls_id"), errors="coerce").fillna(-1).astype(np.int64).to_numpy()
    com["kang_legacy_objectId"] = np.where(ls >= 0, ((ls >> ls_id_shift) & ((1 << (40 - ls_id_shift)) - 1)) * 1_000_000
                                           + (ls & 0xFFFF), -1)
    k = res.merge(com, on="oid", how="left").rename(columns={"oid": "sn_id"})
    k["sn_id"] = k.sn_id.astype(np.int64)
    if k.sn_id.duplicated().any():
        raise SystemExit("duplicate transients in Kang's table")
    return k


def match_positions(cand, sn, ra, dec, rule, band, profile):
    """objectId of the candidate each (ra, dec) position names under ``rule`` (-1: none) and the separation of
    the position from that candidate (arcsec). ``cand`` holds every candidate pair, ``sn`` the transient ids."""
    oid = np.full(len(sn), -1, np.int64)
    sep = np.full(len(sn), np.nan)
    cid = cand.sn_id.to_numpy()
    lo, hi = np.searchsorted(cid, sn, "left"), np.searchsorted(cid, sn, "right")
    if rule == "frac_semimajor":
        a, _, _, q, phi, ok = L.light_model(cand, band, profile)
    cra, cdec, coid = cand.coord_ra.to_numpy(), cand.coord_dec.to_numpy(), cand.objectId.to_numpy()
    for i in range(len(sn)):
        s = slice(lo[i], hi[i])
        if hi[i] == lo[i] or not np.isfinite(ra[i]):
            continue
        x, y = sky_to_plane(cra[s], cdec[s], ra[i], dec[i])
        d = np.hypot(x, y)
        if rule == "frac_semimajor":
            t = np.deg2rad(phi[s])
            u, v = x * np.cos(t) + y * np.sin(t), -x * np.sin(t) + y * np.cos(t)
            r = np.where(ok[s], np.hypot(u, v / q[s]) / a[s], np.inf)
        else:
            r = d
        k = int(np.argmin(r))
        if r[k] <= MATCH[rule]:
            oid[i], sep[i] = coid[s][k], d[k]
    return oid, sep


def build(assoc_path, cand_path, index, rule, cut, workers, ls_id_shift=20):
    """Per-transient table: SBW, DLR, DELIGHT and DELIGHT-SEx hosts (objectId, -1 = hostless) and diagnostics."""
    t = pq.read_table(assoc_path)
    attrs = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if not k.startswith(b"pandas")}
    dp2_lib.install(attrs.get("exclude_suspect", "keep"))
    band, profile = attrs["band"], attrs["profile"]
    assoc = t.to_pandas(ignore_metadata=True)
    kang = read_kang(index, ls_id_shift)
    df = assoc.merge(kang, on="sn_id", how="inner").sort_values("sn_id").reset_index(drop=True)
    cand = L.load_candidates(cand_path, bands=band)
    cand = cand[cand.sn_id.isin(df.sn_id)].reset_index(drop=True)
    legacy = "legacy_type" in pq.read_schema(cand_path).names   # candidates from src/dp2_candidates_legacy.py
    df["kang_objectId"] = df.kang_legacy_objectId if legacy else df.kang_dp2_objectId
    sn = df.sn_id.to_numpy()
    for m, (rc, dc) in {"delight": ("ra_delight", "dec_delight"), "sex": ("ra_sex", "dec_sex")}.items():
        df[f"{m}_objectId"], df[f"{m}_match_sep"] = match_positions(cand, sn, df[rc].to_numpy(float),
                                                                    df[dc].to_numpy(float), rule, band, profile)
        x, y = sky_to_plane(df[rc].to_numpy(float), df[dc].to_numpy(float), df.sn_ra.to_numpy(), df.sn_dec.to_numpy())
        df[f"{m}_sn_sep"] = np.hypot(x, y)
    for r in MATCH:   # every reference_match rule for the sensitivity block of the metrics
        df[f"delight_objectId__{r}"] = match_positions(cand, sn, df.ra_delight.to_numpy(float),
                                                       df.dec_delight.to_numpy(float), r, band, profile)[0]
    # P_SBW, d_DLR and object-table status of the DELIGHT-matched object
    if len(df):
        smp = df[["sn_id", "sn_ra", "sn_dec"]].assign(host_objectId=df.delight_objectId.to_numpy())
        re = associate_all(smp, cand, band=band, profile=profile, search_radius=attrs["search_radius"],
                           dlr_max=attrs["dlr_max"], dlr_radius=attrs["dlr_radius"], extendedness=attrs["extendedness"],
                           quality=attrs["candidate_quality"], workers=workers).set_index("sn_id").loc[sn]
        for c in ("sbw_objectId", "dlr_objectId"):
            bad = int((re[c].to_numpy() != df[c].to_numpy()).sum())
            if bad:
                raise SystemExit(f"recomputed {c} differs from dp2_associations for {bad} transients")
        df["p_sbw_delight"] = re.p_sbw_true.to_numpy()
        df["d_dlr_delight"] = re.d_dlr_true.to_numpy()
        df["delight_status"] = np.where(df.delight_objectId.to_numpy() < 0, "no_candidate",
                                        re.true_host_status.to_numpy())
        smp = df[["sn_id", "sn_ra", "sn_dec"]].assign(host_objectId=df.sex_objectId.to_numpy())
        rx = associate_all(smp, cand, band=band, profile=profile, search_radius=attrs["search_radius"],
                           dlr_max=attrs["dlr_max"], dlr_radius=attrs["dlr_radius"], extendedness=attrs["extendedness"],
                           quality=attrs["candidate_quality"], workers=workers).set_index("sn_id").loc[sn]
        df["p_sbw_sex"] = rx.p_sbw_true.to_numpy()
    else:
        df["p_sbw_delight"], df["d_dlr_delight"], df["delight_status"], df["p_sbw_sex"] = np.nan, np.nan, "", np.nan
    df["sbw_host"] = np.where(sbw_hostless(df, cut), -1, df.sbw_objectId.to_numpy())
    df["dlr_host"] = df.dlr_objectId.to_numpy()
    df["delight_host"], df["sex_host"] = df.delight_objectId.to_numpy(), df.sex_objectId.to_numpy()
    df.attrs = dict(attrs, reference_match=rule, sbw_hostless_cut=cut, catalogue="legacy_dr10" if legacy else "dp2_object")
    return df


def pair_cells(df, a, b, t):
    """Hosted / hostless cells for a pair of methods. Where both assign a host and method a is SBW, 'same' means the
    other method's host is in SBW's credible set (P_SBW >= t) and 'different' that SBW excludes it; for two
    single-pick methods 'same' is the same object."""
    ha, hb = df[f"{a}_host"].to_numpy(), df[f"{b}_host"].to_numpy()
    na, nb = ha < 0, hb < 0
    if a == "sbw":
        same = np.nan_to_num(df[f"p_sbw_{b}"].to_numpy(float), nan=0.0) >= t
    else:
        same = ha == hb
    return {"same": ~na & ~nb & same, "different": ~na & ~nb & ~same, "a_only_hostless": na & ~nb,
            "b_only_hostless": ~na & nb, "both_hostless": na & nb}


def subsets(df):
    yield "all", df
    for s, g in df.groupby("image_survey", sort=True):
        yield s, g


def agreement(df, out, t):
    rows = []
    for name, g in subsets(df):
        for a, b in PAIRS:
            c = pair_cells(g, a, b, t)
            r = {"subset": name, "method_a": LABEL[a], "method_b": LABEL[b], "n": len(g)}
            r.update({k: int(v.sum()) for k, v in c.items()})
            r.update({f"frac_{k}": (float(v.mean()) if len(g) else np.nan) for k, v in c.items()})
            hosted = c["same"] | c["different"]
            r["frac_same_given_both_hosted"] = float(c["same"].sum() / hosted.sum()) if hosted.any() else np.nan
            rows.append(r)
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.6g")


def w(k, n):
    p, lo, hi = wilson(int(k), int(n))
    return {"k": int(k), "n": int(n), "frac": float(p), "lo68": float(lo), "hi68": float(hi)}


def p_bands(g, t):
    """SBW hostless, else P_SBW of the DELIGHT-matched object (0 when the object is not among SBW's
    candidates) in > 0.5 / t-0.5 / < t; over transients whose DELIGHT position matched a candidate."""
    g = g[g.delight_host >= 0]
    sh = g.sbw_host.to_numpy() < 0
    p = np.nan_to_num(g.p_sbw_delight.to_numpy(float), nan=0.0)
    n = len(g)
    return {"n": n, "gt_0p5": w((~sh & (p > 0.5)).sum(), n), "t_to_0p5": w((~sh & (p >= t) & (p <= 0.5)).sum(), n),
            "lt_t": w((~sh & (p < t)).sum(), n), "sbw_hostless": w(sh.sum(), n),
            "mean_p_sbw_delight_sbw_hosted": float(p[~sh].mean()) if (~sh).any() else None,
            "delight_status": {k: int(v) for k, v in g.delight_status.value_counts().items()}}


def metrics(df, out, t):
    radius = L.RADIUS[df.attrs["search_radius"]]
    res = {"settings": {k: df.attrs[k] for k in ("band", "profile", "search_radius", "dlr_max", "dlr_radius",
                                                  "extendedness", "candidate_quality", "reference_match",
                                                  "sbw_hostless_cut", "catalogue")} | {"p_threshold": t},
           "n": len(df), "subsets": {}}
    for name, g in subsets(df):
        n = len(g)
        r = {"n": n, "hostless": {LABEL[m]: w((g[f"{m}_host"] < 0).sum(), n) for m in METHODS}, "pairs": {}}
        for a, b in PAIRS:
            c = pair_cells(g, a, b, t)
            r["pairs"][f"{LABEL[a]}|{LABEL[b]}"] = {k: w(v.sum(), n) for k, v in c.items()} | {
                "same_given_both_hosted": w(c["same"].sum(), (c["same"] | c["different"]).sum()),
                "same_definition": f"P_SBW({LABEL[b]} host) >= {t:g}" if a == "sbw" else "same object"}
        r["p_sbw_of_delight_host"] = p_bands(g, t)
        r["delight_beyond_search_radius"] = w((g.delight_sn_sep > radius).sum(), n)
        r["sex_beyond_search_radius"] = w((g.sex_sn_sep > radius).sum(), n)
        res["subsets"][name] = r
    # the DELIGHT match under every reference_match rule (the decision picks the one used above)
    sh, dh = df.sbw_host.to_numpy(), df.dlr_host.to_numpy()
    res["reference_match_sensitivity"] = {}
    for r in MATCH:
        g = df[f"delight_objectId__{r}"].to_numpy()
        res["reference_match_sensitivity"][r] = {
            "delight_matched": w((g >= 0).sum(), len(df)),
            "sbw_same_as_delight": w(((g >= 0) & (sh == g)).sum(), len(df)),
            "dlr_same_as_delight": w(((g >= 0) & (dh == g)).sum(), len(df)),
            "sbw_same_given_both_hosted": w(((g >= 0) & (sh == g)).sum(), ((g >= 0) & (sh >= 0)).sum())}
    # cross-check of the DELIGHT match against Kang's own nearest DP2 object (within 1" of the CNN position),
    # or his Legacy DR10 counterpart when the candidates are Legacy sources
    k = df.kang_objectId.to_numpy()
    both = (k >= 0) & (df.delight_host.to_numpy() >= 0)
    res["kang_crosscheck"] = {
        "catalogue": df.attrs.get("catalogue", "dp2_object"),
        "kang_matched": w((k >= 0).sum(), len(df)), "ours_matched": w((df.delight_host >= 0).sum(), len(df)),
        "same_object_where_both_matched": w((both & (k == df.delight_host.to_numpy())).sum(), both.sum()),
        "kang_matched_ours_not": int(((k >= 0) & (df.delight_host.to_numpy() < 0)).sum()),
        "ours_matched_kang_not": int(((k < 0) & (df.delight_host.to_numpy() >= 0)).sum()),
        "median_sep_delight_to_matched_arcsec": float(np.nanmedian(df.delight_match_sep)) if len(df) else None}
    res["median_sep_sn_delight_arcsec"] = float(np.nanmedian(df.delight_sn_sep)) if len(df) else None
    res["median_sep_sn_sex_arcsec"] = float(np.nanmedian(df.sex_sn_sep)) if len(df) else None
    with open(out, "w") as f:
        json.dump(res, f, indent=1)


def cases(df, cand_path, out):
    from dp2_cases import ellipses
    sh, dh = df.sbw_host.to_numpy() < 0, df.dlr_host.to_numpy() < 0
    gh, xh = df.delight_host.to_numpy() < 0, df.sex_host.to_numpy() < 0
    s, d, g, x = (df[f"{m}_host"].to_numpy() for m in METHODS)
    classes = {"sbw_ne_delight": ~sh & ~gh & (s != g), "dlr_ne_delight": ~dh & ~gh & (d != g),
               "sbw_hostless_delight_hosted": sh & ~gh, "delight_nocand": gh, "sex_ne_delight": ~xh & ~gh & (x != g)}
    rng = np.random.default_rng(SEED)
    picks = []
    for name, sel in classes.items():
        idx = np.flatnonzero(sel)
        idx = rng.choice(idx, size=min(N_PER_CLASS, len(idx)), replace=False)
        picks.append(df.iloc[np.sort(idx)].assign(**{"class": name, "class_n": int(sel.sum())}))
    c = pd.concat(picks, ignore_index=True)
    band, profile, dlr_radius = df.attrs["band"], df.attrs["profile"], df.attrs["dlr_radius"]
    cand = L.load_candidates(cand_path, bands=band)
    cand = cand[cand.sn_id.isin(c.sn_id)].drop_duplicates(["sn_id", "objectId"]).reset_index(drop=True)
    cand = pd.concat([cand[["sn_id", "objectId", "coord_ra", "coord_dec"]], ellipses(cand, band, profile, dlr_radius)],
                     axis=1).set_index(["sn_id", "objectId"])
    for name, col in [("sbw", "sbw_objectId"), ("dlr", "dlr_objectId")]:
        pos = cand.reindex(pd.MultiIndex.from_arrays([c.sn_id.to_numpy(), c[col].to_numpy()]))
        c[f"{name}_ra"], c[f"{name}_dec"] = pos.coord_ra.to_numpy(), pos.coord_dec.to_numpy()
        for k in ("a", "b", "phi"):
            c[f"{name}_{k}"] = pos[f"{name}_{k}"].to_numpy()
    c["delight_ra"], c["delight_dec"], c["sex_ra"], c["sex_dec"] = c.ra_delight, c.dec_delight, c.ra_sex, c.dec_sex
    c["case_id"] = [f"{k}_{i:02d}" for k, i in zip(c["class"], c.groupby("class").cumcount())]
    c["label"] = [f"{k} | {sid} | {sv} | P(SBW) {p:.2f}, P(DELIGHT host) {pg:.2f} | flb {f:.2g} | "
                  f"DELIGHT {sg:.1f}\" from SN"
                  for k, sid, sv, p, pg, f, sg in zip(c["class"], c.sn_id, c.image_survey, c.p_sbw_max,
                                                      c.p_sbw_delight, c.flb_best, c.delight_sn_sep)]
    cols = (["case_id", "class", "class_n", "sn_id", "sn_ra", "sn_dec", "image_survey", "label"]
            + [f"{n}_{k}" for n in ("sbw", "dlr") for k in ("ra", "dec", "a", "b", "phi")]
            + ["delight_ra", "delight_dec", "sex_ra", "sex_dec", "p_sbw_max", "p_sbw_dlr", "p_sbw_delight",
               "d_dlr_delight", "delight_status", "flb_best", "d_dlr_min", "delight_sn_sep", "sex_sn_sep",
               "delight_match_sep", "std_delight", "hostsize", "sbw_objectId", "dlr_objectId", "delight_objectId",
               "sex_objectId", "kang_dp2_objectId"])
    c[[k for k in cols if k in c.columns]].to_csv(out, index=False, float_format="%.10g")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["agreement", "metrics", "cases"])
    ap.add_argument("--associations", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--delight", required=True, help="index csv (path, size, mtime, sha256) of Kang's two tables")
    ap.add_argument("--reference-match", required=True, choices=sorted(MATCH))
    ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--p-threshold", choices=sorted(P_TRUE), help="agreement and metrics: SBW agrees when P_SBW(other host) >= t")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--kang-ls-id-shift", type=int, default=20,
                    help="bit position of the brick id in Kang's legacy_ls_id (the Kang cross-check on Legacy candidates)")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = build(a.associations, a.candidates, a.delight, a.reference_match, HOSTLESS_CUT[a.hostless_cut], a.workers,
               a.kang_ls_id_shift)
    print(f"{len(df)} transients with DELIGHT predictions", file=sys.stderr)
    if a.mode in ("agreement", "metrics") and not a.p_threshold:
        raise SystemExit("--p-threshold is required for agreement and metrics")
    if a.mode == "agreement":
        agreement(df, a.output, P_TRUE[a.p_threshold])
    elif a.mode == "metrics":
        metrics(df, a.output, P_TRUE[a.p_threshold])
    else:
        cases(df, a.candidates, a.output)


if __name__ == "__main__":
    main()
