"""DP2 (real data) light models, columns and HATS access, plugged into the DP0.2 machinery.

``install()`` replaces the DP0.2 column lists and the light-model / DLR-ellipse functions of ``dp02_lib``
with the DP2 ones, so ``dp02_associate.py``, ``dp02_sweep.py`` and friends run unchanged with ``--tier dp2``
(``associate_one`` looks the functions up in its module at call time, and the worker pool is forked after
the swap).

Conventions
-----------
* DP2 has two multiband MultiProFit fits with one shape for all bands and per-band fluxes: ``exponential_*``
  (n = 1) and ``sersic_*`` (free n). ``{m}_reff_x``, ``{m}_reff_y`` (tract pixels, 0.2") and ``{m}_rho`` are
  the covariance-ellipse parameters in the tract pixel frame (+x West, +y North, the ``sky_to_plane`` frame);
  ``sbw_host.sersic.ellipse_axes`` turns them into semi-axes and an angle counter-clockwise from +x.
  ``{m}_theta`` is never used (its sign convention is odd).
* The fitted model is MultiProFit's Gaussian mixture, so ``reff`` is the mixture scale parameter and the
  real half-light radius is ``half_light_factor(n) * reff`` (1.386 for n = 1); SBW evaluates the mixture
  (``profile="mixture"`` in ``sbw_host.core``), as the ``sbw_host`` package does on DP2.
* sbw_profile on DP2: exponential -> the exponential fit with ``{b}_exponentialFlux``; sersic -> the free-n
  Sersic fit with ``{b}_sersicFlux``; moments -> a Gaussian on the PSF-subtracted HSM moments of the band
  with ``{b}_cModelFlux``, exactly as on DP0.2.
* DLR ``half_light`` is the exponential fit's half-light ellipse whatever the SBW profile, so DLR does not
  move with sbw_profile (on DP0.2 it is the CModel bulge+disk ellipse); ``moments`` is the rms ellipse of the
  band's PSF-subtracted moments, as on DP0.2.
* The ``suspect`` flag is the ``sbw_host`` one (model centroid > 1" from the detection, model / cModel flux
  outside 0.5-2, or brighter than 24 mag yet fainter than 25.5 mag/arcsec^2 inside R_e); with
  exclude_suspect = drop a suspect candidate counts as having no usable light model.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import dp02_lib as L
from sbw_host import associate as SA
from sbw_host import sersic as S
from sbw_host import sersic_mix as SM

PIX = 0.2
FITS = ("exponential", "sersic")
FIT_GENERAL = ("ra", "dec", "reff_x", "reff_y", "rho", "unknown_flag", "no_data_flag", "chi2_reduced")
GENERAL = (["sn_id", "objectId", "coord_ra", "coord_dec", "sep_arcsec", "refExtendedness", "parentObjectId",
            "parent_nChild", "tract", "patch", "detect_fromBlend", "detect_isIsolated", "deblend_blendNChild",
            "shape_xx", "shape_yy", "shape_xy", "shape_flag", "sersic_index", "exponential_delta_lnL_fit_ps"]
           + [f"{m}_{k}" for m in FITS for k in FIT_GENERAL])
PER_BAND = ["exponentialFlux", "exponentialFluxErr", "sersicFlux", "sersicFluxErr", "cModelFlux", "cModelFluxErr",
            "cModel_flag", "blendedness", "ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF", "extendedness"]
# Object-table names that differ from the candidate-table names
RENAME = {"deblend_parentNChild": "parent_nChild"}
PROFILE_LABEL = {"exponential": "exponential fit (n = 1, MultiProFit mixture)",
                 "sersic": "free-n Sersic fit (MultiProFit mixture)",
                 "moments": "Gaussian (PSF-subtracted moments)"}
CFG = {"exclude_suspect": "keep"}

_moment_light_model = L.light_model   # the DP0.2 Gaussian-on-moments branch is reused as is


def object_columns(bands=L.BANDS):
    """Object-table columns the candidate extraction reads."""
    inv = {v: k for k, v in RENAME.items()}
    gen = [inv.get(c, c) for c in GENERAL if c not in ("sn_id", "sep_arcsec")]
    return gen + [f"{b}_{c}" for b in bands for c in PER_BAND]


def load_candidates(path, bands=L.BANDS):
    cols = GENERAL + [f"{b}_{c}" for b in bands for c in PER_BAND]
    df = pq.read_table(path, columns=cols).to_pandas(ignore_metadata=True)
    return df.sort_values(["sn_id", "sep_arcsec"], kind="stable").reset_index(drop=True)


def fit_ellipse(c, fit):
    """Half-light semi-major axis (arcsec), axis ratio, angle (deg from +x toward +y), Sersic index and a
    validity mask of one MultiProFit fit."""
    g = lambda k: c[f"{fit}_{k}"].to_numpy(dtype=float)
    rx, ry, rho = g("reff_x") * PIX, g("reff_y") * PIX, g("rho")
    n = c["sersic_index"].to_numpy(dtype=float) if fit == "sersic" else np.ones(len(c))
    ok = np.isfinite(rx) & np.isfinite(ry) & np.isfinite(rho) & (rx > 0) & (ry > 0) & np.isfinite(n) & (n > 0)
    for f in ("unknown_flag", "no_data_flag"):
        ok &= ~c[f"{fit}_{f}"].fillna(True).to_numpy(dtype=bool)
    major, minor, phi = S.ellipse_axes(np.where(ok, rx, 1.0), np.where(ok, ry, 1.0), np.where(ok, rho, 0.0))
    hl = SM.half_light_factor(np.where(ok, n, 1.0))
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.clip(minor / major, 1e-3, 1.0)
    return major * hl, q, phi, n, ok


def suspect(c, band, fit):
    """sbw_host's suspect flag for the ``fit`` model in ``band`` (see the module docstring)."""
    g = lambda k: c[k].to_numpy(dtype=float)
    du, dv = S.sky_offsets(g(f"{fit}_ra"), g(f"{fit}_dec"), g("coord_ra"), g("coord_dec"))
    shift = np.where(np.isfinite(du), np.hypot(du, dv), 0.0)
    flux = g(f"{band}_{fit}Flux")
    a, q, _, _, _ = fit_ellipse(c, fit)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = flux / g(f"{band}_cModelFlux")
        mag = -2.5 * np.log10(flux) + 31.4
        sb = mag + 0.753 + 2.5 * np.log10(np.pi * a * a * q)
    lo, hi = SA.SUSPECT_FLUX_RATIO_RANGE
    return ((shift > SA.SUSPECT_MAX_CENTROID_SHIFT_ARCSEC) | (np.isfinite(ratio) & ((ratio < lo) | (ratio > hi)))
            | ((sb > SA.SUSPECT_DIFFUSE_SB_MAG_ARCSEC2) & (mag < SA.SUSPECT_DIFFUSE_MAX_MAG)))


def light_model(c, band, profile):
    """DP2 version of ``dp02_lib.light_model``: semi-major half-light radius (arcsec), flux (nJy), n, q, phi and
    a validity mask."""
    if profile == "moments":
        return _moment_light_model(c, band, profile)
    if profile not in FITS:
        raise ValueError(f"profile {profile!r} is not one of {L.PROFILES} (stamps need images)")
    a, q, phi, n, ok = fit_ellipse(c, profile)
    flux = c[f"{band}_{profile}Flux"].to_numpy(dtype=float)
    ok &= np.isfinite(flux) & (flux > 0)
    if CFG["exclude_suspect"] == "drop":
        ok &= ~suspect(c, band, profile)
    return a, flux, n, q, phi, ok


def dlr_ellipse(c, band, dlr_radius, q, l1, l2):
    """DLR semi-axes (arcsec), validity and the ellipse angle (deg from +x toward +y)."""
    if dlr_radius == "half_light":
        a, qe, phi, _, ok = fit_ellipse(c, "exponential")
        return a, a * qe, ok, phi
    if dlr_radius == "moments":
        g = lambda k: c[f"{band}_{k}"].to_numpy(dtype=float)
        _, phi, m1, m2, valid = L.moment_shape(g("ixx"), g("iyy"), g("ixy"), g("ixxPSF"), g("iyyPSF"), g("ixyPSF"))
        a, b = np.sqrt(m1) * PIX, np.sqrt(m2) * PIX
        return a, b, valid & np.isfinite(a) & (a > 0) & (b > 0), phi
    raise ValueError(dlr_radius)


def install(exclude_suspect="keep"):
    """Make dp02_lib associate DP2 candidates."""
    CFG["exclude_suspect"] = exclude_suspect
    L.PER_BAND, L.GENERAL = PER_BAND, GENERAL
    L.light_model, L.dlr_ellipse, L.load_candidates = light_model, dlr_ellipse, load_candidates
    L.PROFILE_LABEL.update(PROFILE_LABEL)
    L.SBW_CURVE.update({"exponential": "mixture", "sersic": "mixture"})


# ---------------------------------------------------------------- HATS object_collection
def hats_root(index_csv):
    """Catalogue directory (the one holding ``dataset/``) from the dp2_object index of partition files."""
    p = pd.read_csv(index_csv, nrows=1)["path"].iloc[0]
    return p[: p.index("/dataset/")]


def partition_boxes(root):
    """Per-partition path and coord_ra / coord_dec bounding box from the HATS per-partition statistics."""
    st = pd.read_parquet(os.path.join(root, "per_partition_statistics.parquet"),
                         columns=["Norder", "Npix", "coord_ra::min_value", "coord_ra::max_value",
                                  "coord_dec::min_value", "coord_dec::max_value"])
    st.columns = ["Norder", "Npix", "ra_min", "ra_max", "dec_min", "dec_max"]
    st["path"] = [os.path.join(root, "dataset", f"Norder={o}", f"Dir={(p // 10000) * 10000}", f"Npix={p}.parquet")
                  for o, p in zip(st.Norder, st.Npix)]
    return st


def partitions_near(boxes, ra, dec, radius_arcsec):
    """Indices into ``boxes`` of the partitions whose bounding box, padded by the radius, holds (ra, dec);
    a box spanning > 180 deg in RA (one wrapping through 0) is always kept."""
    pad = radius_arcsec / 3600.0
    padra = pad / max(np.cos(np.deg2rad(min(abs(dec) + pad, 89.9))), 1e-3)
    wrap = (boxes.ra_max - boxes.ra_min).to_numpy() > 180
    ok = ((boxes.dec_min.to_numpy() - pad <= dec) & (dec <= boxes.dec_max.to_numpy() + pad)
          & (wrap | ((boxes.ra_min.to_numpy() - padra <= ra) & (ra <= boxes.ra_max.to_numpy() + padra))))
    return np.flatnonzero(ok)


def unit_vectors(ra, dec):
    r, d = np.deg2rad(np.asarray(ra, float)), np.deg2rad(np.asarray(dec, float))
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def chord(arcsec):
    return 2 * np.sin(np.deg2rad(arcsec / 3600.0) / 2)
