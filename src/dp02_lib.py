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
PROFILES = ("exponential", "sersic", "moments")
PROFILE_LABEL = {"exponential": "exponential (n = 1, CModel disk radius)",
                 "sersic": "de Vaucouleurs (n = 4, CModel bulge radius)",
                 "moments": "Gaussian (PSF-subtracted moments)"}
SERSIC_N = {"exponential": 1.0, "sersic": 4.0, "moments": 0.5}
COLORS = {"sbw": "#2a6fdb", "dlr": "#e07b39", "agree": "#555555", "ref": "#999999"}

PER_BAND = ["cModelFlux", "bdReB", "bdReD", "bdFluxB", "bdFluxD",
            "ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF", "extendedness"]
GENERAL = ["sn_id", "objectId", "coord_ra", "coord_dec", "sep_arcsec", "refExtendedness"]


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
    """Semi-major half-light radius (arcsec), flux (nJy), n, q, phi and a validity mask."""
    g = lambda k: c[f"{band}_{k}"].to_numpy(dtype=float)
    flux = g("cModelFlux")
    q, phi, l1, l2, valid = moment_shape(g("ixx"), g("iyy"), g("ixy"), g("ixxPSF"), g("iyyPSF"), g("ixyPSF"))
    if profile == "exponential":
        a = g("bdReD") * PIX / np.sqrt(q)
    elif profile == "sersic":
        a = g("bdReB") * PIX / np.sqrt(q)
    elif profile == "moments":
        a = GAUSS_HALF_LIGHT * np.sqrt(l1) * PIX
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


def associate_one(sn_ra, sn_dec, host_objectId, c, band, profile, radius, dlr_max, dlr_radius, extendedness):
    """Associate one SN against its candidate rows ``c`` (a DataFrame slice)."""
    out = {"n_cand": len(c), "n_cand_used": 0, "true_in_cand": False,
           "sbw_objectId": -1, "p_sbw_max": np.nan, "p_sbw_true": np.nan, "rank_true": -1,
           "flb_best": np.nan, "flb_true": np.nan, "m_true": np.nan, "sep_sbw": np.nan,
           "dlr_objectId": -1, "d_dlr_min": np.nan, "d_dlr_true": np.nan, "d_dlr_sbw": np.nan, "sep_dlr": np.nan,
           "sep_true": np.nan}
    if len(c) == 0:
        return out
    keep = c["sep_arcsec"].to_numpy() <= radius
    if extendedness == "ext05":
        ext = c["refExtendedness"].to_numpy(dtype=float)
        keep &= ext >= 0.5
    c = c[keep]
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
    # ---- SBW
    if ok.any():
        res = sbw_model(0.0, 0.0, x[ok], y[ok], r_eff=a[ok], flux=flux[ok], n=n, q=q[ok], phi=phi[ok])
        p = np.full(len(c), 0.0); p[ok] = np.nan_to_num(res.p_sbw)
        flb = np.full(len(c), np.nan); flb[ok] = res.frac_light_beyond
        m = np.full(len(c), np.nan); m[ok] = res.m
        k = int(np.argmax(p))
        if p[k] > 0:
            out["sbw_objectId"] = int(oid[k]); out["p_sbw_max"] = p[k]; out["flb_best"] = flb[k]
            out["sep_sbw"] = sep[k]
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
                  dlr_radius, extendedness, workers=None) -> pd.DataFrame:
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
                     dlr_radius=dlr_radius, extendedness=extendedness)
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
                     dlr_radius=dlr_radius, extendedness=extendedness).items():
        out.attrs[k] = v
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
