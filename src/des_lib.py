"""DES DR2 (SExtractor) light models and DLR ellipses, plugged into the DP0.2 machinery.

``install()`` replaces the column lists and the light-model / DLR-ellipse functions of ``dp02_lib`` with the DES
ones (as ``dp2_lib.install`` does for DP2), so ``dp02_associate.py``, ``dp02_sweep.py`` and ``dp2_fig.py`` run
unchanged through ``src/des_run.py``, which also restricts their band loops to DES's griz(y).

Conventions
-----------
* Shapes are the SExtractor detection-image (r+i+z) moments ``a_image``, ``b_image`` (rms semi-axes, 0.263"
  pixels) and ``theta_j2000`` (deg East of North), one shape for every band, as DES-SN used them. In the
  sky_to_plane frame (+x West, +y North) the major axis lies at ``phi = theta_j2000 + 90`` deg counter-clockwise
  from +x (checked on DES DR2 images, output des_angle_check). des_candidates writes them as the DP0.2 moment
  columns ``{b}_ixx/iyy/ixy`` in arcsec^2 with zero PSF moments, so ``dp02_lib.moment_shape`` returns
  q = B/A, phi and l1 = A^2, l2 = B^2 (arcsec^2) unchanged; they are not PSF-subtracted.
* The half-light size is the band's ``flux_radius`` (circular radius enclosing half the AUTO flux, PSF
  included), taken as the determinant radius sqrt(a b) of the half-light ellipse, so the semi-major axis is
  flux_radius / sqrt(q), as for DP0.2's bdReD. Fluxes are ``flux_auto`` (mag = 30 - 2.5 log10 flux) in nJy.
* sbw_profile on DES (DR2 has no Sersic fit): exponential -> n = 1 on the flux_radius ellipse; sersic -> n = 4
  (de Vaucouleurs) on the same ellipse, as DP0.2's n = 4 with the bulge radius; moments -> a Gaussian whose
  rms ellipse is A, B (half-light radius 1.177 A).
* DLR ``moments`` is the A, B, THETA ellipse (the DES-SN definition); ``half_light`` the flux_radius ellipse
  of the band.
* Point sources: refExtendedness = 1 for extended_class_coadd >= 2 ("mostly galaxies"), 0 for 0 / 1 (stars,
  candidate stars), NaN (cut by ext05) for -9 (no measurement).
* candidate_quality reads cModelFlux / cModelFluxErr = flux_auto / fluxerr_auto, cModel_flag = SExtractor
  flags >= 4 or imaflags_iso > 0, blendedness NaN (no such column: strict = basic + the sibling cut, and DES DR2
  has no parent ids, so the sibling cut never fires).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import dp02_lib as L

PIX = 0.263                     # arcsec per DES coadd pixel
NJY = 10 ** (0.4 * (31.4 - 30))  # flux_auto (zero point 30) -> nJy
BANDS = "grizy"
GENERAL = ["sn_id", "objectId", "coord_ra", "coord_dec", "sep_arcsec", "refExtendedness", "parentObjectId",
           "parent_nChild", "a_image", "b_image", "theta_j2000", "kron_radius", "extended_class_coadd"]
PER_BAND = ["cModelFlux", "cModelFluxErr", "cModel_flag", "blendedness", "ixx", "iyy", "ixy", "ixxPSF", "iyyPSF",
            "ixyPSF", "extendedness", "flux_radius", "mag_auto", "class_star", "spread_model", "spreaderr_model",
            "flags", "imaflags_iso"]
PROFILE_LABEL = {"exponential": "exponential (n = 1, FLUX_RADIUS ellipse)",
                 "sersic": "de Vaucouleurs (n = 4, FLUX_RADIUS ellipse)",
                 "moments": "Gaussian (SExtractor A, B)"}
SERSIC_N = {"exponential": 1.0, "sersic": 4.0, "moments": 0.5}


def load_candidates(path, bands=BANDS):
    cols = GENERAL + [f"{b}_{c}" for b in bands for c in PER_BAND]
    df = pq.read_table(path, columns=cols).to_pandas(ignore_metadata=True)
    return df.sort_values(["sn_id", "sep_arcsec"], kind="stable").reset_index(drop=True)


def shape(c):
    """rms semi-axes A, B (arcsec), axis ratio, angle (deg from +x West toward +y North) and validity."""
    A = c["a_image"].to_numpy(dtype=float) * PIX
    B = c["b_image"].to_numpy(dtype=float) * PIX
    ok = np.isfinite(A) & np.isfinite(B) & (A > 0) & (B > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(ok, np.clip(B / A, 0.05, 1.0), 1.0)
    phi = c["theta_j2000"].to_numpy(dtype=float) + 90.0
    return A, B, q, np.where(np.isfinite(phi), phi, 0.0), ok & np.isfinite(phi)


def light_model(c, band, profile):
    """DES version of ``dp02_lib.light_model``: semi-major half-light radius (arcsec), flux (nJy), n, q, phi and
    a validity mask."""
    A, _, q, phi, ok = shape(c)
    flux = c[f"{band}_cModelFlux"].to_numpy(dtype=float)
    if profile in ("exponential", "sersic"):
        a = c[f"{band}_flux_radius"].to_numpy(dtype=float) * PIX / np.sqrt(q)
    elif profile == "moments":
        a = L.GAUSS_HALF_LIGHT * A
    else:
        raise ValueError(f"profile {profile!r} is not one of {L.PROFILES} (stamps need images)")
    ok = ok & np.isfinite(a) & (a > 0) & np.isfinite(flux) & (flux > 0)
    return a, flux, SERSIC_N[profile], q, phi, ok


def dlr_ellipse(c, band, dlr_radius, q, l1, l2):
    """DLR semi-axes (arcsec), validity and angle (deg from +x toward +y)."""
    A, B, qs, phi, ok = shape(c)
    if dlr_radius == "moments":
        return A, B, ok, phi
    if dlr_radius == "half_light":
        a = c[f"{band}_flux_radius"].to_numpy(dtype=float) * PIX / np.sqrt(qs)
        return a, a * qs, ok & np.isfinite(a) & (a > 0), phi
    raise ValueError(dlr_radius)


def install():
    """Make dp02_lib associate DES DR2 candidates."""
    L.PER_BAND, L.GENERAL, L.BANDS = PER_BAND, GENERAL, BANDS
    L.light_model, L.dlr_ellipse, L.load_candidates = light_model, dlr_ellipse, load_candidates
    L.PROFILE_LABEL.update(PROFILE_LABEL)
    L.SERSIC_N.update(SERSIC_N)


# ---------------------------------------------------------------- catalogue -> candidate schema
def candidate_columns(obj: pd.DataFrame) -> pd.DataFrame:
    """Rewrite des_dr2.main rows (one per object) into the candidate schema above, without sn_id / sep_arcsec."""
    out = pd.DataFrame({"objectId": obj.coadd_object_id.astype(np.int64),
                        "coord_ra": obj.ra.astype(float), "coord_dec": obj.dec.astype(float)})
    ecc = obj.extended_class_coadd.to_numpy(dtype=float)
    out["refExtendedness"] = np.where(ecc < 0, np.nan, (ecc >= 2).astype(float))
    out["parentObjectId"] = np.int64(0)
    out["parent_nChild"] = np.int64(0)
    for k in ("a_image", "b_image", "theta_j2000", "kron_radius", "extended_class_coadd"):
        out[k] = obj[k].to_numpy()
    A, B, _, phi, ok = shape(obj)
    t = np.deg2rad(phi)
    ixx = np.where(ok, A ** 2 * np.cos(t) ** 2 + B ** 2 * np.sin(t) ** 2, np.nan)
    iyy = np.where(ok, A ** 2 * np.sin(t) ** 2 + B ** 2 * np.cos(t) ** 2, np.nan)
    ixy = np.where(ok, (A ** 2 - B ** 2) * np.sin(t) * np.cos(t), np.nan)
    for b in BANDS:
        f = obj[f"flux_auto_{b}"].to_numpy(dtype=float)
        out[f"{b}_cModelFlux"] = f * NJY
        out[f"{b}_cModelFluxErr"] = obj[f"fluxerr_auto_{b}"].to_numpy(dtype=float) * NJY
        out[f"{b}_cModel_flag"] = (obj[f"flags_{b}"].to_numpy() >= 4) | (obj[f"imaflags_iso_{b}"].to_numpy() > 0)
        out[f"{b}_blendedness"] = np.nan
        out[f"{b}_ixx"], out[f"{b}_iyy"], out[f"{b}_ixy"] = ixx, iyy, ixy
        for k in ("ixxPSF", "iyyPSF", "ixyPSF"):
            out[f"{b}_{k}"] = 0.0
        out[f"{b}_extendedness"] = out["refExtendedness"]
        out[f"{b}_flux_radius"] = obj[f"flux_radius_{b}"].to_numpy(dtype=float)
        for k in ("mag_auto", "class_star", "spread_model", "spreaderr_model", "flags", "imaflags_iso"):
            out[f"{b}_{k}"] = obj[f"{k}_{b}"].to_numpy()
    return out


def unit_vectors(ra, dec):
    r, d = np.deg2rad(np.asarray(ra, float)), np.deg2rad(np.asarray(dec, float))
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def chord(arcsec):
    return 2 * np.sin(np.deg2rad(np.asarray(arcsec, float) / 3600.0) / 2)


def arcsec_from_chord(ch):
    return np.rad2deg(2 * np.arcsin(np.clip(ch / 2, 0, 1))) * 3600.0
