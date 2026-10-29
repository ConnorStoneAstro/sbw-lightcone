"""Legacy Surveys DR10 (south) Tractor sources as DP2 candidates, written in the dp2_candidates schema.

Called by ``src/dp2_candidates.py --catalogue legacy_dr10``. Every sweep tile of the legacy_dr10_sweeps
index is matched against all positions (a position needing a tile that is not indexed is an error), so a
transient near a tile edge gets the sources of both tiles. TYPE DUP rows (Gaia sources fitted as another
source, zero flux) are dropped.

The rows are rewritten so that ``dp2_lib`` reads them as if they were the DP2 MultiProFit fits, and
``dp02_associate.py --tier dp2`` and everything after it run unchanged:

* objectId = BRICKID * 1_000_000 + OBJID; tract = BRICKID, patch = -1; brickname, legacy_type, maskbits,
  nobs_r, fracflux_r and legacy_sersic are kept as extra columns.
* refExtendedness (and every ``{b}_extendedness``) = 0 for TYPE PSF, 1 otherwise, so extendedness_cut ext05
  removes PSF sources and ``none`` keeps them.
* Shape: SHAPE_R is the half-light radius along the *semi-major* axis (tractor ``EllipseE.getRaDecBasis``:
  major axis r_e, minor r_e b/a), q = (1 - |e|) / (1 + |e|), and the position angle East of North is
  PA = 0.5 atan2(SHAPE_E2, SHAPE_E1) (checked on r-band cutouts, output dp2_legacy_angle_check). In the DP2
  frame (+x West, +y North, angle from +x toward +y) that is phi = PA + 90 deg.
* Profiles: ``exponential_*`` is the Tractor ellipse with n = 1 for every source; ``sersic_*`` the same
  ellipse with the Tractor index (REX, EXP 1; DEV 4; SER the fitted SERSIC; PSF 0.5). SBW evaluates the
  MultiProFit Gaussian-mixture Sersic of that index, so reff = R_e / half_light_factor(n) is written as the
  covariance-ellipse parameters reff_x, reff_y, rho in 0.2" pixels. PSF sources (no shape) get a round
  Gaussian with R_e = PSFSIZE_R / 2 (the half-light radius of a Gaussian is its HWHM) and are flagged
  no-data when NOBS_R = 0.
* Fluxes: FLUX_{g,r,i,z} (nanomaggies) * 3631 = nJy as exponentialFlux, sersicFlux and cModelFlux (Tractor
  has one model flux per band); errors from FLUX_IVAR; u and y are NaN. blendedness = FRACFLUX of the band
  (flux fraction from neighbours, the nearest Tractor analogue). There are no HSM moments, so the moments
  light model and dlr_radius moments are unusable on this catalogue (all NaN).
"""
from __future__ import annotations

import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.spatial import cKDTree

import dp2_lib as D
from legacy_index_sweeps import tiles_near
from sbw_host import sersic_mix as SM

NMGY_NJY = 3631.0
BAND_COL = {"g": "G", "r": "R", "i": "I", "z": "Z"}
TYPE_N = {"REX": 1.0, "EXP": 1.0, "DEV": 4.0, "PSF": 0.5}
READ = ["BRICKID", "BRICKNAME", "OBJID", "TYPE", "RA", "DEC", "DCHISQ", "SHAPE_R", "SHAPE_E1", "SHAPE_E2",
        "SERSIC", "PSFSIZE_R", "NOBS_R", "RCHISQ_R", "MASKBITS"] + [
    f"{k}_{B}" for B in BAND_COL.values() for k in ("FLUX", "FLUX_IVAR", "FRACFLUX")]
_G = {}


def covariance_params(a, q, phi_deg):
    """reff_x, reff_y, rho of the ellipse with semi-major a, axis ratio q and angle phi (from +x to +y);
    the inverse of sbw_host.sersic.ellipse_axes."""
    b, c, s = a * q, np.cos(np.deg2rad(phi_deg)), np.sin(np.deg2rad(phi_deg))
    sx, sy = np.sqrt(a * a * c * c + b * b * s * s), np.sqrt(a * a * s * s + b * b * c * c)
    return sx, sy, (a * a - b * b) * c * s / (sx * sy)


def _match(path):
    sn = _G["sn"]
    with fits.open(path, memmap=True) as h:
        t = h[1].data
        tree = cKDTree(D.unit_vectors(t["RA"], t["DEC"]))
        hits = tree.query_ball_point(D.unit_vectors(sn.sn_ra, sn.sn_dec), D.chord(_G["radius"]))
        pairs = [(i, j) for i, js in enumerate(hits) for j in js]
        if not pairs:
            return None
        si, oj = np.array(pairs).T
        rows = np.unique(oj)
        native = lambda x: np.asarray(x).astype(np.asarray(x).dtype.newbyteorder("="))   # FITS is big-endian
        cat = pd.DataFrame({k: (native(t[k][rows]) if k != "DCHISQ" else list(native(t[k][rows]))) for k in READ})
    cat = cat.iloc[np.searchsorted(rows, oj)].reset_index(drop=True)
    cat.insert(0, "sn_id", sn.sn_id.to_numpy()[si])
    cat.insert(1, "sn_ra", sn.sn_ra.to_numpy()[si])
    cat.insert(2, "sn_dec", sn.sn_dec.to_numpy()[si])
    return cat


def to_dp2_schema(cat):
    """Rewrite matched sweep rows as dp2_candidates rows (see the module docstring)."""
    typ = cat.TYPE.astype(str).str.strip().to_numpy()
    cat = cat[typ != "DUP"].reset_index(drop=True)
    typ = typ[typ != "DUP"]
    f = lambda k: cat[k].to_numpy(dtype=float)
    psf = typ == "PSF"
    e = np.hypot(f("SHAPE_E1"), f("SHAPE_E2"))
    q = np.where(psf | (typ == "REX"), 1.0, np.clip((1 - e) / (1 + e), 1e-3, 1.0))
    pa = 0.5 * np.degrees(np.arctan2(f("SHAPE_E2"), f("SHAPE_E1")))
    phi = np.where(q < 1, pa + 90.0, 0.0)
    re = np.where(psf, f("PSFSIZE_R") / 2.0, f("SHAPE_R"))
    n = np.array([TYPE_N.get(t, np.nan) for t in typ])
    n = np.where(typ == "SER", f("SERSIC"), n)
    nodata = ~(np.isfinite(re) & (re > 0)) | (psf & (f("NOBS_R") <= 0))
    re = np.where(nodata, np.nan, re)
    a = {"sn_id": cat.sn_id.to_numpy(np.int64),
         "objectId": cat.BRICKID.to_numpy(np.int64) * 1_000_000 + cat.OBJID.to_numpy(np.int64),
         "coord_ra": f("RA"), "coord_dec": f("DEC")}
    u, v = D.unit_vectors(a["coord_ra"], a["coord_dec"]), D.unit_vectors(cat.sn_ra, cat.sn_dec)
    a["sep_arcsec"] = np.degrees(2 * np.arcsin(np.clip(np.linalg.norm(u - v, axis=1) / 2, 0, 1))) * 3600
    ext = np.where(psf, 0.0, 1.0)
    a["refExtendedness"] = ext
    a["parentObjectId"], a["parent_nChild"] = np.int64(0), 0
    a["tract"], a["patch"] = cat.BRICKID.to_numpy(np.int64), -1
    a["detect_fromBlend"], a["detect_isIsolated"], a["deblend_blendNChild"] = False, np.nan, 0
    a["shape_xx"] = a["shape_yy"] = a["shape_xy"] = np.nan
    a["shape_flag"] = True
    a["sersic_index"] = n
    dchisq = np.vstack(cat.DCHISQ.to_numpy()).astype(float)   # PSF, REX, DEV, EXP, SER
    a["exponential_delta_lnL_fit_ps"] = 0.5 * (dchisq[:, 3] - dchisq[:, 0])
    for fit, nn in (("exponential", np.ones(len(cat))), ("sersic", n)):
        hl = SM.half_light_factor(np.where(np.isfinite(nn) & (nn > 0), nn, 1.0))
        rx, ry, rho = covariance_params(re / hl / D.PIX, q, phi)
        a[f"{fit}_ra"], a[f"{fit}_dec"] = a["coord_ra"], a["coord_dec"]
        a[f"{fit}_reff_x"], a[f"{fit}_reff_y"], a[f"{fit}_rho"] = rx, ry, rho
        a[f"{fit}_unknown_flag"] = ~np.isfinite(nn)
        a[f"{fit}_no_data_flag"] = nodata
        a[f"{fit}_chi2_reduced"] = f("RCHISQ_R")
    for b in "ugrizy":
        B = BAND_COL.get(b)
        flux = f(f"FLUX_{B}") * NMGY_NJY if B else np.full(len(cat), np.nan)
        with np.errstate(divide="ignore"):
            err = np.where(f(f"FLUX_IVAR_{B}") > 0, NMGY_NJY / np.sqrt(f(f"FLUX_IVAR_{B}")), np.nan) if B else flux
        for k in ("exponentialFlux", "sersicFlux", "cModelFlux"):
            a[f"{b}_{k}"], a[f"{b}_{k}Err"] = flux, err
        a[f"{b}_cModel_flag"] = ~np.isfinite(err)
        a[f"{b}_blendedness"] = f(f"FRACFLUX_{B}") if B else np.nan
        for k in ("ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF"):
            a[f"{b}_{k}"] = np.nan
        a[f"{b}_extendedness"] = ext
    a["brickname"] = cat.BRICKNAME.astype(str).str.strip().to_numpy()
    a["legacy_type"], a["legacy_sersic"] = typ, f("SERSIC")
    a["maskbits"], a["nobs_r"], a["fracflux_r"] = cat.MASKBITS.to_numpy(np.int64), cat.NOBS_R.to_numpy(np.int64), f("FRACFLUX_R")
    return pd.DataFrame(a)


def build(sn, sweeps_csv, radius, workers):
    """dp2_candidates rows for positions ``sn`` (sn_id, sn_ra, sn_dec) from the indexed sweep tiles."""
    paths = pd.read_csv(sweeps_csv)["path"].tolist()
    have = {os.path.basename(p) for p in paths}
    need = set().union(*[tiles_near(r, d, radius / 3600.0) for r, d in zip(sn.sn_ra, sn.sn_dec)])
    if need - have:
        raise SystemExit(f"sweep tiles {sorted(need - have)} are needed but not in {sweeps_csv}")
    paths = [p for p in paths if os.path.basename(p) in need]
    print(f"{len(sn)} positions, {len(paths)} sweep tiles to read", file=sys.stderr)
    _G.update(sn=sn, radius=radius)
    with Pool(min(workers, len(paths))) as pool:
        parts = [p for p in pool.map(_match, paths) if p is not None]   # ordered: deterministic output
    return to_dp2_schema(pd.concat(parts, ignore_index=True))
