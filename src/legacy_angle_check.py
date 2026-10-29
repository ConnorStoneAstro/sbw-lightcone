"""Check the Legacy Surveys DR10 Tractor ellipse convention against the images.

Tractor stores a galaxy shape as SHAPE_R (half-light radius, arcsec) and SHAPE_E1, SHAPE_E2 with
|e| = (1 - q) / (1 + q). The candidate builder (src/dp2_candidates_legacy.py) needs the position angle;
this script measures it on Legacy DR10 r-band cutouts (viewer cutout.fits, ls-dr10, 0.262"/pix) of large,
elongated, unmasked galaxies from one sweep tile and compares it with theta = 0.5 atan2(e2, e1) under the
two sign conventions. The image angle comes from Gaussian-weighted second moments (weight sigma = SHAPE_R,
background = median of the cutout border) transformed to (East, North) with the cutout's CD matrix, so it
does not depend on how the cutout is oriented. Page 1: residual histograms PA_image - (+theta) and
PA_image - (-theta), with <cos 2 delta>; page 2: cutouts with the catalogue ellipse drawn for each sign.
"""
from __future__ import annotations

import argparse
import io
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

URL = "https://www.legacysurvey.org/viewer/cutout.fits"
PIX = 0.262


def select(path, n, seed):
    with fits.open(path, memmap=True) as h:
        t = h[1].data
        typ = np.char.strip(t["TYPE"].astype(str))
        r, e1, e2 = t["SHAPE_R"], t["SHAPE_E1"], t["SHAPE_E2"]
        e = np.hypot(e1, e2)
        ok = (np.isin(typ, ["EXP", "DEV", "SER"]) & (r > 5) & (r < 15) & (e > 0.35) & (t["MASKBITS"] == 0)
              & (t["FLUX_R"] > 100) & (t["NOBS_R"] > 0))
        idx = np.flatnonzero(ok)
        df = pd.DataFrame({k: np.asarray(t[k][idx]).astype(float) for k in
                           ("RA", "DEC", "SHAPE_R", "SHAPE_E1", "SHAPE_E2", "FLUX_R")})
        df["TYPE"] = typ[idx]
    print(f"{len(df)} large elongated galaxies in {path}", file=sys.stderr)
    return df.sample(min(n, len(df)), random_state=seed).reset_index(drop=True)


def cutout(ra, dec, size, session):
    for k in range(4):
        r = session.get(URL, params={"ra": ra, "dec": dec, "layer": "ls-dr10", "pixscale": PIX, "size": size,
                                     "bands": "r"}, timeout=60)
        if r.status_code == 200:
            h = fits.open(io.BytesIO(r.content))
            data = np.squeeze(h[0].data).astype(float)
            return data, WCS(h[0].header).celestial
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"HTTP {r.status_code}")


def image_pa(data, wcs, ra, dec, re_arcsec):
    """Position angle (deg East of North) and axis ratio of the Gaussian-weighted second moments."""
    ny, nx = data.shape
    x0, y0 = wcs.world_to_pixel_values(ra, dec)
    border = np.concatenate([data[:3].ravel(), data[-3:].ravel(), data[:, :3].ravel(), data[:, -3:].ravel()])
    d = data - np.nanmedian(border)
    yy, xx = np.mgrid[:ny, :nx]
    dx, dy = xx - x0, yy - y0
    w = np.exp(-0.5 * (dx * dx + dy * dy) / (re_arcsec / PIX) ** 2) * np.nan_to_num(d)
    m = np.array([[np.sum(w * dx * dx), np.sum(w * dx * dy)], [np.sum(w * dx * dy), np.sum(w * dy * dy)]]) / w.sum()
    # Jacobian (East, North) per pixel, measured from the WCS itself at the galaxy
    (r0, d0), (r1, d1), (r2, d2) = np.array(
        wcs.pixel_to_world_values(np.array([x0, x0 + 1, x0]), np.array([y0, y0, y0 + 1]))).T
    cosd = np.cos(np.deg2rad(d0))
    J = np.array([[(r1 - r0) * cosd, (r2 - r0) * cosd], [d1 - d0, d2 - d0]])
    s = J @ m @ J.T   # (East, North) moments
    th_en = 0.5 * np.degrees(np.arctan2(2 * s[0, 1], s[0, 0] - s[1, 1]))   # major axis from East toward North
    ev = np.linalg.eigvalsh(s)
    return (90.0 - th_en) % 180.0, float(np.sqrt(ev[0] / ev[1]))


def wrap90(x):
    return (np.asarray(x) + 90.0) % 180.0 - 90.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweeps", required=True, help="legacy_dr10_sweeps index csv")
    ap.add_argument("--tile", default="sweep-050m030-055m025.fits")
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--min-interval", type=float, default=1.0)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    paths = pd.read_csv(a.sweeps)["path"]
    path = paths[paths.str.endswith(a.tile)].iloc[0]
    g = select(path, a.n, a.seed)
    theta = 0.5 * np.degrees(np.arctan2(g.SHAPE_E2, g.SHAPE_E1))
    e = np.hypot(g.SHAPE_E1, g.SHAPE_E2)
    g["q_cat"] = (1 - e) / (1 + e)
    s = requests.Session()
    pa, qi, imgs = [], [], []
    for _, row in g.iterrows():
        size = int(np.clip(8 * row.SHAPE_R / PIX, 64, 400))
        try:
            data, wcs = cutout(row.RA, row.DEC, size, s)
            p, q = image_pa(data, wcs, row.RA, row.DEC, row.SHAPE_R)
        except Exception as ex:   # noqa: BLE001
            print(f"cutout {row.RA:.5f} {row.DEC:.5f}: {ex}", file=sys.stderr)
            data, wcs, p, q = None, None, np.nan, np.nan
        pa.append(p), qi.append(q), imgs.append((data, wcs))
        time.sleep(a.min_interval)
    g["pa_image"], g["q_image"], g["theta"] = pa, qi, theta
    ok = np.isfinite(g.pa_image) & (g.q_image < 0.85)   # a measurable image elongation
    res = {"+theta": wrap90(g.pa_image - g.theta), "-theta": wrap90(g.pa_image + g.theta),
           "90-theta": wrap90(g.pa_image - (90 - g.theta)), "90+theta": wrap90(g.pa_image - (90 + g.theta))}
    summ = {k: (float(np.mean(np.cos(np.deg2rad(2 * v[ok])))), float(np.median(np.abs(v[ok])))) for k, v in res.items()}
    for k, (c2, med) in summ.items():
        print(f"PA_image vs PA = {k}: <cos 2 delta> {c2:+.3f}, median |delta| {med:.1f} deg (n {ok.sum()})",
              file=sys.stderr)
    with PdfPages(a.output) as pdf:
        fig, axs = plt.subplots(1, 4, figsize=(14, 3.6), sharey=True)
        for ax, (k, v) in zip(axs, res.items()):
            ax.hist(v[ok], bins=np.arange(-90, 91, 7.5), color="0.4")
            ax.set_title(f"PA(image) - PA, PA = {k}\n<cos 2d> {summ[k][0]:+.2f}, median |d| {summ[k][1]:.0f} deg",
                         fontsize=9)
            ax.set_xlabel("delta (deg)")
        axs[0].set_ylabel(f"galaxies (n {ok.sum()}; q_image < 0.85)")
        fig.suptitle("Legacy DR10 Tractor ellipse angle vs r-band cutout moments; "
                     "theta = 0.5 atan2(SHAPE_E2, SHAPE_E1) in deg, PA East of North", fontsize=10)
        fig.tight_layout()
        pdf.savefig(fig), plt.close(fig)
        fig, axs = plt.subplots(3, 4, figsize=(12, 9.5))
        for ax, (_, row), (data, wcs) in zip(axs.ravel(), g.iterrows(), imgs):
            if data is None:
                ax.axis("off")
                continue
            ny, nx = data.shape
            yy, xx = np.mgrid[-0.5:ny, -0.5:nx]
            cra, cdec = wcs.wcs_pix2world(xx, yy, 0)
            ex = (cra - row.RA) * np.cos(np.deg2rad(row.DEC)) * 3600
            no = (cdec - row.DEC) * 3600
            v = np.nan_to_num(data)
            lo, hi = np.percentile(v, [5, 99.7])
            ax.pcolormesh(ex, no, np.arcsinh(np.clip((v - lo) / (hi - lo), 0, 1) * 20), cmap="gray_r",
                          shading="flat", rasterized=True)
            for sign, col in ((+1, "C0"), (-1, "C3")):
                pa_deg = sign * row.theta   # PA East of North; matplotlib angle from +x (East here) = 90 - PA
                ax.add_patch(Ellipse((0, 0), 4 * row.SHAPE_R, 4 * row.SHAPE_R * row.q_cat, angle=90 - pa_deg,
                                     fill=False, ec=col, lw=1.2))
            h = 4 * row.SHAPE_R
            ax.set_xlim(h, -h), ax.set_ylim(-h, h), ax.set_aspect("equal")
            ax.set_xticks([]), ax.set_yticks([])
            ax.set_title(f"{row.TYPE} R_e {row.SHAPE_R:.1f}\" q {row.q_cat:.2f} | img {row.pa_image:.0f}, "
                         f"theta {row.theta:.0f}", fontsize=7)
        fig.suptitle("blue: PA = +theta, red: PA = -theta (2 R_e ellipse); East left, North up", fontsize=10)
        fig.tight_layout()
        pdf.savefig(fig), plt.close(fig)


if __name__ == "__main__":
    main()
