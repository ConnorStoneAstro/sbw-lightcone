"""Shared DP0.2 machinery: candidate loading, light models from DP0.2 columns,
per-SN SBW + DLR association, and a parallel batch runner.

Conventions
-----------
* DP0.2 ``bdReB`` / ``bdReD`` are determinant radii ``(xx*yy - xy^2)^(1/4)`` of
  the CModel de Vaucouleurs / exponential half-light ellipse in pixels (0.2"),
  i.e. ``sqrt(a*b)``; the semi-major axis is ``r / sqrt(q)``.
* The v23 ``bdE1`` column is broken (missing parentheses in the functor), so
  axis ratio and orientation come from the PSF-subtracted HSM source moments
  ``ixx - ixxPSF`` etc., in the tract pixel frame (+x West, +y North), which is
  the frame ``sbw_host.core.sky_to_plane`` uses.
* Fluxes are nJy; intensities are nJy / arcsec^2.
"""
from __future__ import annotations

import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from sbw_host.core import sbw_model, sky_to_plane
from sbw_host.dlr import d_dlr

PIX = 0.2                      # arcsec per DP0.2 coadd pixel
GAUSS_HALF_LIGHT = 1.17741     # half-light radius of a Gaussian in units of sigma
BANDS = "ugrizy"
RADIUS = {"r15": 15.0, "r30": 30.0, "r60": 60.0}
DLR_MAX = {"d4": 4.0, "d5": 5.0}
HOSTLESS_CUT = {"flb001": 0.001, "flb01": 0.01, "flb02": 0.02, "flb05": 0.05}   # SBW hostless when flb_best < cut
PROFILES = ("exponential", "sersic", "moments")
PROFILE_LABEL = {"exponential": "exponential (n = 1, CModel disk radius)",
                 "sersic": "de Vaucouleurs (n = 4, CModel bulge radius)",
                 "moments": "Gaussian (PSF-subtracted moments)"}
SERSIC_N = {"exponential": 1.0, "sersic": 4.0, "moments": 0.5}
COLORS = {"sbw": "#2a6fdb", "dlr": "#e07b39", "agree": "#555555", "ref": "#999999"}

PER_BAND = ["cModelFlux", "bdReB", "bdReD", "bdFluxB", "bdFluxD",
            "ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF", "extendedness",
            "cModelFluxErr", "cModel_flag", "blendedness"]
GENERAL = ["sn_id", "objectId", "coord_ra", "coord_dec", "sep_arcsec", "refExtendedness",
           "parentObjectId", "parent_nChild"]

# candidate_quality: object-table cuts applied to candidates before association (headline band)
QUALITY = ("none", "basic", "strict")
SNR_MIN = 5.0          # basic: cModelFlux / cModelFluxErr
BLEND_MAX = 0.5        # strict: blendedness above this is dropped
# strict: more siblings (same parent among the SN's candidates within the search radius) than this is dropped.
# Sibling counts barely separate true hosts from wrong picks on DP0.2 (baseline r30, 27 Oct 2026: 96% of true
# hosts are deblended children, median 10 siblings vs 12 for wrong SBW picks; > 30 removes 8.7% of true hosts and
# 14.5% of wrong picks); 50 is the smallest value that removes < 1% of true hosts (0.85%, wrong picks 2.4%).
SIBLING_MAX = 50
P_TRUE = {"p001": 0.01, "p005": 0.05}   # below this P_SBW the true host is effectively excluded
P_TOP_MIN = 1e-3       # p_top keeps every candidate probability above this
# true_host_status: why the true host is or is not among the candidates SBW and DLR see, first reason wins
STATUS = ("absent", "point_source", "unusable", "quality", "kept")


def load_candidates(path, bands=BANDS):
    cols = GENERAL + [f"{b}_{c}" for b in bands for c in PER_BAND]
    t = pq.read_table(path, columns=cols)
    df = t.to_pandas(ignore_metadata=True)
    return df.sort_values(["sn_id", "sep_arcsec"], kind="stable").reset_index(drop=True)


def moment_shape(ixx, iyy, ixy, pxx, pyy, pxy):
    """Axis ratio, position angle (deg, from +x toward +y) and principal second
    moments (pixel^2) of the PSF-subtracted moments; falls back to the raw
    moments where subtraction is not positive definite."""
    def eig(mxx, myy, mxy):
        tr = mxx + myy
        d = np.sqrt((mxx - myy) ** 2 + 4 * mxy ** 2)
        l1, l2 = (tr + d) / 2, (tr - d) / 2
        ok = np.isfinite(l1) & np.isfinite(l2) & (l2 > 0) & (l1 > 0)
        return l1, l2, ok
    l1, l2, ok = eig(ixx - pxx, iyy - pyy, ixy - pxy)
    r1, r2, rok = eig(ixx, iyy, ixy)
    use_raw = ~ok & rok
    l1 = np.where(use_raw, r1, l1)
    l2 = np.where(use_raw, r2, l2)
    mxx = np.where(use_raw, ixx, ixx - pxx)
    myy = np.where(use_raw, iyy, iyy - pyy)
    mxy = np.where(use_raw, ixy, ixy - pxy)
    valid = ok | rok
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.sqrt(np.clip(l2 / l1, 0.0025, 1.0))
        phi = 0.5 * np.degrees(np.arctan2(2 * mxy, mxx - myy))
    q = np.where(valid, q, 1.0)
    phi = np.where(valid, phi, 0.0)
    return q, phi, l1, l2, valid


def light_model(c, band, profile):
    """Semi-major half-light radius (arcsec), flux (nJy), n, q, phi and a validity mask. A candidate
    without a usable radius or flux for the profile is an object-table failure and is dropped upstream."""
    g = lambda k: c[f"{band}_{k}"].to_numpy(dtype=float)
    flux = g("cModelFlux")
    q, phi, l1, l2, valid = moment_shape(g("ixx"), g("iyy"), g("ixy"), g("ixxPSF"), g("iyyPSF"), g("ixyPSF"))
    a_mom = GAUSS_HALF_LIGHT * np.sqrt(l1) * PIX
    if profile == "exponential":
        a = g("bdReD") * PIX / np.sqrt(q)
    elif profile == "sersic":
        a = g("bdReB") * PIX / np.sqrt(q)
    elif profile == "moments":
        a = a_mom
    else:
        raise ValueError(f"profile {profile!r} is not one of {PROFILES} (stamps need images, not available here)")
    ok = valid & np.isfinite(a) & (a > 0) & np.isfinite(flux) & (flux > 0)
    return a, flux, SERSIC_N[profile], q, phi, ok


def dlr_ellipse(c, band, dlr_radius, q, l1, l2):
    """Semi-axes (arcsec) of the DLR ellipse."""
    g = lambda k: c[f"{band}_{k}"].to_numpy(dtype=float)
    if dlr_radius == "half_light":
        fb, fd = np.clip(g("bdFluxB"), 0, None), np.clip(g("bdFluxD"), 0, None)
        rb, rd = g("bdReB"), g("bdReD")
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (fb * rb + fd * rd) / (fb + fd)
        r = np.where(np.isfinite(r), r, np.where(np.isfinite(rd), rd, rb))
        a = r * PIX / np.sqrt(q)
        b = a * q
    elif dlr_radius == "moments":
        a, b = np.sqrt(l1) * PIX, np.sqrt(l2) * PIX
    else:
        raise ValueError(dlr_radius)
    ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    return a, b, ok


def siblings(parent):
    """Number of other rows in ``parent`` with the same non-zero parentObjectId (0 = not deblended)."""
    parent = np.nan_to_num(np.asarray(parent, dtype=float), nan=0.0).astype(np.int64)
    _, inv, cnt = np.unique(parent, return_inverse=True, return_counts=True)
    return np.where(parent > 0, cnt[inv] - 1, 0)


def quality_mask(c, band, level, n_sib):
    """True for candidates that pass the candidate_quality cut, from object-table columns only."""
    keep = np.ones(len(c), bool)
    if level == "none":
        return keep
    g = lambda k: c[f"{band}_{k}"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = g("cModelFlux") / g("cModelFluxErr")
    keep &= ~c[f"{band}_cModel_flag"].fillna(True).to_numpy(dtype=bool) & (snr >= SNR_MIN)
    if level == "strict":
        if SIBLING_MAX is None:
            raise ValueError("SIBLING_MAX is not set in dp02_lib.py")
        keep &= ~(g("blendedness") > BLEND_MAX) & (n_sib <= SIBLING_MAX)
    elif level != "basic":
        raise ValueError(level)
    return keep


def associate_one(sn_ra, sn_dec, host_objectId, c, band, profile, radius, dlr_max, dlr_radius, extendedness,
                  quality="none"):
    """Associate one SN against its candidate rows ``c`` (a DataFrame slice).

    The candidate list is cut in order: search radius, extendedness_cut (point sources), no usable radius or
    flux for the profile (always), candidate_quality; ``true_host_status`` records which step removed the
    true host (``absent`` when it is not a primary object within the search radius) and n_cut_* how many
    candidates each step removed."""
    out = {"n_cand": len(c), "n_cand_radius": 0, "n_cut_point_source": 0, "n_cut_unusable": 0, "n_cut_quality": 0,
           "n_cand_used": 0, "true_in_cand": False, "true_host_status": "absent", "n_sib_true": -1,
           "sbw_objectId": -1, "p_sbw_max": np.nan, "p_sbw_true": np.nan, "rank_true": -1, "n_sib_sbw": -1,
           "p_top": [],
           "flb_best": np.nan, "flb_max": np.nan, "flb_true": np.nan, "m_true": np.nan, "sep_sbw": np.nan,
           "dlr_objectId": -1, "d_dlr_min": np.nan, "d_dlr_true": np.nan, "d_dlr_sbw": np.nan, "sep_dlr": np.nan,
           "sep_true": np.nan}
    if len(c) == 0:
        return out
    c = c[c["sep_arcsec"].to_numpy() <= radius]
    out["n_cand_radius"] = len(c)
    if len(c) == 0:
        return out
    n_sib = siblings(c["parentObjectId"].to_numpy())
    istrue = c["objectId"].to_numpy() == host_objectId
    if istrue.any():
        out["true_host_status"] = "kept"
        out["n_sib_true"] = int(n_sib[istrue][0])
    steps = []
    if extendedness == "ext05":
        steps.append(("point_source", c["refExtendedness"].to_numpy(dtype=float) >= 0.5))
    ok_all = light_model(c, band, profile)[-1]
    steps.append(("unusable", ok_all))
    steps.append(("quality", quality_mask(c, band, quality, n_sib)))
    keep = np.ones(len(c), bool)
    for name, m in steps:
        drop = keep & ~m
        out[f"n_cut_{name}"] = int(drop.sum())
        if (drop & istrue).any():
            out["true_host_status"] = name
        keep &= m
    c, n_sib = c[keep], n_sib[keep]
    oid = c["objectId"].to_numpy()
    out["n_cand_used"] = len(c)
    if len(c) == 0:
        return out
    it = np.flatnonzero(oid == host_objectId)
    it = int(it[0]) if len(it) else -1
    out["true_in_cand"] = it >= 0
    sep = c["sep_arcsec"].to_numpy(dtype=float)
    if it >= 0:
        out["sep_true"] = sep[it]
    x, y = sky_to_plane(c["coord_ra"].to_numpy(), c["coord_dec"].to_numpy(), sn_ra, sn_dec)
    a, flux, n, q, phi, ok = light_model(c, band, profile)
    _, _, l1, l2, _ = moment_shape(*(c[f"{band}_{k}"].to_numpy(dtype=float)
                                     for k in ("ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF")))
    # ---- SBW (every remaining candidate has a usable light model)
    res = sbw_model(0.0, 0.0, x, y, r_eff=a, flux=flux, n=n, q=q, phi=phi)
    p = np.nan_to_num(np.asarray(res.p_sbw, dtype=float))
    flb = np.asarray(res.frac_light_beyond, dtype=float)
    m = np.asarray(res.m, dtype=float)
    out["flb_max"] = float(np.nanmax(flb)) if np.isfinite(flb).any() else np.nan
    out["p_top"] = np.sort(p[p > P_TOP_MIN])[::-1].astype(np.float32).tolist()
    k = int(np.argmax(p))
    if p[k] > 0:
        out["sbw_objectId"] = int(oid[k]); out["p_sbw_max"] = p[k]; out["flb_best"] = flb[k]
        out["sep_sbw"] = sep[k]; out["n_sib_sbw"] = int(n_sib[k])
    if it >= 0:
        out["p_sbw_true"] = p[it]; out["flb_true"] = flb[it]; out["m_true"] = m[it]
        out["rank_true"] = int((p > p[it]).sum()) + 1
    # ---- DLR
    da, db, dok = dlr_ellipse(c, band, dlr_radius, q, l1, l2)
    d = np.full(len(c), np.inf)
    d[dok] = d_dlr(-x[dok], -y[dok], da[dok], db[dok], phi[dok])
    if np.isfinite(d).any():
        k = int(np.argmin(d))
        out["d_dlr_min"] = d[k]
        if d[k] <= dlr_max:
            out["dlr_objectId"] = int(oid[k]); out["sep_dlr"] = sep[k]
    if it >= 0 and np.isfinite(d[it]):
        out["d_dlr_true"] = d[it]
    if out["sbw_objectId"] >= 0:
        ks = np.flatnonzero(oid == out["sbw_objectId"])[0]
        out["d_dlr_sbw"] = d[ks] if np.isfinite(d[ks]) else np.nan
    return out


_G = {}


def associate_all(sample: pd.DataFrame, cand: pd.DataFrame, *, band, profile, search_radius, dlr_max,
                  dlr_radius, extendedness, quality="none", workers=None) -> pd.DataFrame:
    """Associate every SN of ``sample`` using its rows in ``cand`` (both keyed by sn_id)."""
    workers = workers or int(os.environ.get("DP02_WORKERS", "16"))
    t0 = time.time()
    sample = sample.sort_values("sn_id").reset_index(drop=True)
    cand = cand.sort_values(["sn_id", "sep_arcsec"], kind="stable").reset_index(drop=True)
    cid = cand["sn_id"].to_numpy()
    sid = sample["sn_id"].to_numpy()
    starts = np.searchsorted(cid, sid, side="left")
    ends = np.searchsorted(cid, sid, side="right")
    # rows for SN i are cand[starts[i]:ends[i]] (empty when the SN has no candidates)
    _G["sample"] = {k: sample[k].to_numpy() for k in ("sn_id", "sn_ra", "sn_dec", "host_objectId")}
    _G["cand"] = cand
    _G["starts"] = np.column_stack([starts, ends])
    _G["cfg"] = dict(band=band, profile=profile, radius=RADIUS[search_radius], dlr_max=DLR_MAX[dlr_max],
                     dlr_radius=dlr_radius, extendedness=extendedness, quality=quality)
    n = len(sample)
    chunks = [(lo, min(lo + 2000, n)) for lo in range(0, n, 2000)]
    rows = []
    with Pool(workers) as pool:
        for part in pool.imap(_run_chunk_pairs, chunks, chunksize=1):
            rows.extend(part)
    out = pd.DataFrame(rows)
    out = sample.merge(out, on="sn_id", how="left")
    out["sbw_correct"] = (out["sbw_objectId"] >= 0) & (out["sbw_objectId"] == out["host_objectId"])
    out["dlr_correct"] = (out["dlr_objectId"] >= 0) & (out["dlr_objectId"] == out["host_objectId"])
    out["dlr_hostless"] = out["dlr_objectId"] < 0
    out["agree"] = (out["sbw_objectId"] >= 0) & (out["sbw_objectId"] == out["dlr_objectId"])
    for k, v in dict(band=band, profile=profile, search_radius=search_radius, dlr_max=dlr_max,
                     dlr_radius=dlr_radius, extendedness=extendedness, candidate_quality=quality).items():
        out.attrs[k] = v
    # what the object-table cuts removed: candidates (pairs within the search radius) and true hosts
    nr = max(int(out.n_cand_radius.sum()), 1)
    for k in ("point_source", "unusable", "quality"):
        out.attrs[f"frac_candidates_cut_{k}"] = float(out[f"n_cut_{k}"].sum() / nr)
    hosted = out.host_objectId.to_numpy() >= 0
    if hosted.any():
        st = out.true_host_status.to_numpy()[hosted]
        for k in STATUS:
            out.attrs[f"frac_true_hosts_{k}"] = float((st == k).mean())
    print(f"associated {n} SNe ({band}, {profile}) in {time.time() - t0:.0f}s", file=sys.stderr)
    return out


def _run_chunk_pairs(idx_range):
    lo, hi = idx_range
    S, C, se, cfg = _G["sample"], _G["cand"], _G["starts"], _G["cfg"]
    rows = []
    for i in range(lo, hi):
        c = C.iloc[se[i, 0]:se[i, 1]]
        r = associate_one(S["sn_ra"][i], S["sn_dec"][i], S["host_objectId"][i], c, **cfg)
        r["sn_id"] = S["sn_id"][i]
        rows.append(r)
    return rows


def sbw_hostless(df, cut, col="flb_best"):
    """SBW declares the transient hostless: no SBW host, or its frac_light_beyond below ``cut``."""
    f = df[col].to_numpy(dtype=float)
    return (df["sbw_objectId"].to_numpy() < 0) | ~(f >= cut)


def dlr_hostless(df, t):
    """DLR declares the transient hostless: no candidate with d_DLR <= t."""
    return ~(df["d_dlr_min"].to_numpy(dtype=float) <= t)


def matrix_cells(df, cut, truth_col="host_objectId"):
    """3x3 SBW (true / wrong / hostless at flb < cut) x DLR (true / wrong / hostless) counts against
    ``truth_col``, plus how many of the both-wrong SNe have SBW and DLR on the same object."""
    t = df[truth_col].to_numpy()
    s, d = df["sbw_objectId"].to_numpy(), df["dlr_objectId"].to_numpy()
    hl = sbw_hostless(df, cut)
    rows = {"sbw_true": (s == t) & (t >= 0) & ~hl, "sbw_wrong": ~((s == t) & (t >= 0)) & ~hl, "sbw_hostless": hl}
    cols = {"dlr_true": (d == t) & (t >= 0), "dlr_wrong": (d >= 0) & (d != t), "dlr_hostless": d < 0}
    out = {f"{r}/{c}": int((rm & cm).sum()) for r, rm in rows.items() for c, cm in cols.items()}
    out["both_wrong_same_object"] = int((rows["sbw_wrong"] & cols["dlr_wrong"] & (s == d)).sum())
    out["n"] = len(df)
    return out


def p_true(df):
    """P_SBW of the true host, 0 where it is not among the surviving candidates."""
    return np.nan_to_num(df["p_sbw_true"].to_numpy(dtype=float), nan=0.0)


def prob_rows(df, cut, t):
    """SBW outcome rows of the probabilistic matrix: hostless (flb < cut) first, otherwise the band of
    P_SBW(true host): > 0.5, t..0.5, < t (excluded, including a true host that is not a candidate)."""
    hl = sbw_hostless(df, cut)
    p = p_true(df)
    return {"p_gt_0p5": ~hl & (p > 0.5), "p_t_to_0p5": ~hl & (p >= t) & (p <= 0.5),
            "p_lt_t": ~hl & (p < t), "sbw_hostless": hl}


def dlr_cols(df, truth_col="host_objectId"):
    t, d = df[truth_col].to_numpy(), df["dlr_objectId"].to_numpy()
    return {"dlr_true": (d == t) & (t >= 0), "dlr_wrong": (d >= 0) & (d != t), "dlr_hostless": d < 0}


def prob_matrix_cells(df, cut, t):
    """4x3 counts: SBW P_SBW(true) band / hostless x DLR true / wrong / hostless."""
    rows, cols = prob_rows(df, cut, t), dlr_cols(df)
    out = {f"{r}/{c}": int((rm & cm).sum()) for r, rm in rows.items() for c, cm in cols.items()}
    out.update(n=len(df), p_true_threshold=t, sbw_hostless_cut=cut)
    return out


MISS = ("success", "genuine", "cut", "absent")
MISS_LABEL = {"success": "success", "genuine": "cat3: genuine association failure",
              "cut": "cat2: true host removed by object-table cuts",
              "absent": "cat1: true host absent from the catalogue"}


def miss_categories(df, cut, t):
    """Per-SN outcome for SBW and DLR (one of MISS) plus the sub-category, from an association table.

    absent: the true host is not a primary object within the search radius; cut: removed by the extendedness
    cut, by having no usable radius / flux, or by candidate_quality; genuine: the true host survived but SBW
    gives it P < t or declares the SN hostless (flb < cut), resp. DLR assigns another candidate or nothing."""
    st = df["true_host_status"].to_numpy()
    absent, kept = st == "absent", st == "kept"
    cutm = ~absent & ~kept
    hl = sbw_hostless(df, cut)
    p = p_true(df)
    sbw = np.where(absent, "absent", np.where(cutm, "cut", np.where(hl | (p < t), "genuine", "success")))
    sbw_sub = np.where(absent, "absent", np.where(cutm, st, np.where(hl, "sbw_hostless",
                                                                      np.where(p < t, "p_below_threshold", "success"))))
    d, tr = df["dlr_objectId"].to_numpy(), df["host_objectId"].to_numpy()
    dlr_unusable = kept & ~np.isfinite(df["d_dlr_true"].to_numpy(dtype=float))
    dlr = np.where(absent, "absent", np.where(cutm | dlr_unusable, "cut",
                                              np.where((d == tr) & (tr >= 0), "success", "genuine")))
    dlr_sub = np.where(absent, "absent", np.where(cutm, st, np.where(dlr_unusable, "unusable",
                       np.where(dlr == "success", "success", np.where(d >= 0, "other_closer", "beyond_dlr_max")))))
    return {"sbw": sbw, "sbw_sub": sbw_sub, "dlr": dlr, "dlr_sub": dlr_sub}


def random_cells(rnd, cut):
    """2x2 SBW (hostless at flb < cut = correct / assigns a host) x DLR (hostless / assigns) counts for
    random positions, which have no host."""
    hl = sbw_hostless(rnd, cut)
    dh = rnd["dlr_objectId"].to_numpy() < 0
    rows = {"sbw_hostless": hl, "sbw_host": ~hl}
    cols = {"dlr_hostless": dh, "dlr_host": ~dh}
    out = {f"{r}/{c}": int((rm & cm).sum()) for r, rm in rows.items() for c, cm in cols.items()}
    out.update(n=len(rnd), sbw_rejected=float(hl.mean()), dlr_rejected=float(dh.mean()))
    return out


# ---------------------------------------------------------------- cosmoDC2 host truth shapes
BULGE_DOMINATED_BT = 0.5   # bulge_to_total_ratio_i above which a host counts as bulge-dominated


def host_offsets(df, pa_convention="east_of_north"):
    """SN offset from the truth host centre in units of the host's half-light semi-major axis, measured
    along the truth ellipse (elliptical radius), for the bulge, the disk and the dominant component.

    ``df`` needs sn_ra, sn_dec, host_ra, host_dec and the cosmoDC2 columns size_(minor_)bulge_true,
    size_(minor_)disk_true (arcsec), position_angle_true (deg) and bulge_to_total_ratio_i.

    "east_of_north" reads position_angle_true as the major axis East of North, which is how DC2 placed
    the SNe (offset angles align with it, <cos 2(psi - PA)> = 0.59 for q < 0.3 hosts). The *rendered*
    galaxies are mirrored: DP0.2 HSM moments of large elongated hosts give PA = -position_angle_true
    (median |dPA| 0.8 deg), and SN angles are uncorrelated with that axis. SN radii follow the semi-major
    R_e exactly (median r/R_e = 1.00) on a rounder ellipse than the galaxy (axis ratio ~ q^0.4)."""
    x, y = sky_to_plane(df["sn_ra"].to_numpy(), df["sn_dec"].to_numpy(),
                        df["host_ra"].to_numpy(), df["host_dec"].to_numpy())   # arcsec, +x West, +y North
    e, n = -x, y
    th = np.deg2rad(df["position_angle_true"].to_numpy(dtype=float))
    if pa_convention == "east_of_north":
        u, v = e * np.sin(th) + n * np.cos(th), e * np.cos(th) - n * np.sin(th)
    else:   # north_of_east
        u, v = e * np.cos(th) + n * np.sin(th), -e * np.sin(th) + n * np.cos(th)
    out = {}
    for comp in ("bulge", "disk"):
        a = df[f"size_{comp}_true"].to_numpy(dtype=float)
        b = df[f"size_minor_{comp}_true"].to_numpy(dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            q = np.where(a > 0, b / a, np.nan)
            out[comp] = np.where(a > 0, np.hypot(u, v / q) / a, np.nan)
    bd = df["bulge_to_total_ratio_i"].to_numpy(dtype=float) > BULGE_DOMINATED_BT
    out["own"] = np.where(bd, out["bulge"], out["disk"])
    out["bulge_dominated"] = bd
    return out


# ---------------------------------------------------------------- statistics helpers
def wilson(k, n, z=1.0):
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k / n
        den = 1 + z ** 2 / n
        centre = (p + z ** 2 / (2 * n)) / den
        half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return p, centre - half, centre + half


CAL_EDGES = np.array([0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 0.995, 1.0000001])


def calibration_curve(df, edges=CAL_EDGES):
    """Bin SNe by p_sbw_max; return mean p, fraction correct with Wilson band, and counts."""
    d = df[np.isfinite(df["p_sbw_max"])]
    p = d["p_sbw_max"].to_numpy()
    c = d["sbw_correct"].to_numpy()
    idx = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, len(edges) - 2)
    n = np.bincount(idx, minlength=len(edges) - 1)
    k = np.bincount(idx, weights=c, minlength=len(edges) - 1)
    pm = np.bincount(idx, weights=p, minlength=len(edges) - 1) / np.maximum(n, 1)
    frac, lo, hi = wilson(k, n)
    ece = float(np.nansum(np.abs(frac - pm) * n) / max(n.sum(), 1))
    return pd.DataFrame({"p_mean": pm, "frac_correct": frac, "lo": lo, "hi": hi, "n": n}), ece
