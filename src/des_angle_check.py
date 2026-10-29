"""Check the DES DR2 SExtractor ellipse convention against the images.

src/des_lib.py reads ``theta_j2000`` as the major-axis position angle East of North, i.e. phi = theta_j2000 + 90
deg counter-clockwise from +x = West in the sky_to_plane frame. This script measures the angle on DES DR2
r-band coadd cutouts (Data Lab SIA, as src/des_gallery.py) of large, elongated, unflagged galaxies from the
des_dr2_objects cache and compares it with +theta_j2000 and with the mirrored -theta_j2000. The image angle comes
from Gaussian-weighted second moments (weight sigma = the A_IMAGE rms radius, background = median of the cutout
border) transformed to (East, North) with the cutout's own WCS, so it does not depend on how the cutout is
oriented. Page 1: residual histograms with <cos 2 delta> and the median |delta|; page 2: cutouts with the
catalogue ellipse drawn for both signs.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
import pyvo  # noqa: E402
import requests  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

import des_gallery as DG  # noqa: E402
import dp02_gallery as G  # noqa: E402

PIX = 0.263
SEED = 20261029


def select(path, n):
    cols = ["coadd_object_id", "ra", "dec", "a_image", "b_image", "theta_j2000", "extended_class_coadd", "flux_auto_r",
            "flags_r", "imaflags_iso_r"]
    o = pq.read_table(path, columns=cols).to_pandas()
    q = o.b_image / o.a_image
    ok = ((o.extended_class_coadd == 3) & (o.a_image * PIX > 1.5) & (o.a_image * PIX < 6) & (q < 0.5)
          & (o.flags_r == 0) & (o.imaflags_iso_r == 0) & (o.flux_auto_r > 3000))
    print(f"{ok.sum()} large elongated galaxies in the cache", file=sys.stderr)
    return o[ok].sample(min(n, int(ok.sum())), random_state=SEED).reset_index(drop=True)


def image_pa(data, wcs, ra, dec, sigma_arcsec):
    """Position angle (deg East of North) and axis ratio of the Gaussian-weighted second moments."""
    ny, nx = data.shape
    x0, y0 = wcs.world_to_pixel_values(ra, dec)
    border = np.concatenate([data[:3].ravel(), data[-3:].ravel(), data[:, :3].ravel(), data[:, -3:].ravel()])
    d = data - np.nanmedian(border)
    yy, xx = np.mgrid[:ny, :nx]
    dx, dy = xx - x0, yy - y0
    w = np.exp(-0.5 * (dx * dx + dy * dy) / (2 * sigma_arcsec / PIX) ** 2) * np.nan_to_num(d)
    m = np.array([[np.sum(w * dx * dx), np.sum(w * dx * dy)], [np.sum(w * dx * dy), np.sum(w * dy * dy)]]) / w.sum()
    (r0, d0), (r1, d1), (r2, d2) = np.array(
        wcs.pixel_to_world_values(np.array([x0, x0 + 1, x0]), np.array([y0, y0, y0 + 1]))).T
    cosd = np.cos(np.deg2rad(d0))
    J = np.array([[(r1 - r0) * cosd, (r2 - r0) * cosd], [d1 - d0, d2 - d0]])
    s = J @ m @ J.T   # (East, North) moments
    th_en = 0.5 * np.degrees(np.arctan2(2 * s[0, 1], s[0, 0] - s[1, 1]))   # from East toward North
    ev = np.linalg.eigvalsh(s)
    return (90.0 - th_en + 180.0) % 180.0, float(np.sqrt(max(ev[0], 0) / ev[1]))


def wrap(d):
    return (np.asarray(d) + 90.0) % 180.0 - 90.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--objects", required=True)
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--size", type=float, default=30.0)
    ap.add_argument("--cache-dir", default="des_cutouts")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    gal = select(a.objects, a.n)
    sia, session, limiter = pyvo.dal.SIAService(DG.SIA), requests.Session(), G.RateLimiter(0.7)
    rows, imgs = [], []
    for g in gal.itertuples():
        try:
            data, wcs = G.read_image(DG.fetch(g.ra, g.dec, "r", a.size, sia, session, Path(a.cache_dir), limiter))
        except Exception as e:  # noqa: BLE001
            print(f"{g.coadd_object_id}: {e}", file=sys.stderr)
            continue
        pa_img, q_img = image_pa(data, wcs, g.ra, g.dec, g.a_image * PIX)
        rows.append({"objectId": g.coadd_object_id, "pa_image": pa_img, "q_image": q_img, "theta_j2000": g.theta_j2000,
                     "q_cat": g.b_image / g.a_image, "a_arcsec": g.a_image * PIX})
        imgs.append((g, data, wcs))
    r = pd.DataFrame(rows)
    dp, dm = wrap(r.pa_image - r.theta_j2000), wrap(r.pa_image + r.theta_j2000)
    stat = lambda d: (float(np.mean(np.cos(np.deg2rad(2 * d)))), float(np.median(np.abs(d))))
    (cp, mp), (cm, mm) = stat(dp), stat(dm)
    verdict = ("theta_j2000 is East of North (phi = theta_j2000 + 90 from West)" if cp > 0.9 and cp > cm else
               "theta_j2000 is mirrored (PA = -theta_j2000)" if cm > 0.9 else "inconclusive")
    print(f"N {len(r)}: +theta <cos2d> {cp:.3f} median|d| {mp:.1f}; -theta {cm:.3f} {mm:.1f}; {verdict}", file=sys.stderr)
    with PdfPages(a.output) as pdf:
        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        bins = np.linspace(-90, 90, 37)
        ax.hist(dp, bins=bins, histtype="step", lw=1.5, color="#2a6fdb",
                label=f"PA$_{{img}}$ - theta_j2000: <cos 2$\\delta$> = {cp:.3f}, median |$\\delta$| = {mp:.1f} deg")
        ax.hist(dm, bins=bins, histtype="step", lw=1.5, color="#e07b39",
                label=f"PA$_{{img}}$ + theta_j2000: <cos 2$\\delta$> = {cm:.3f}, median |$\\delta$| = {mm:.1f} deg")
        ax.set_xlabel("residual [deg]"); ax.set_ylabel("galaxies"); ax.legend(fontsize=7, loc="upper left")
        ax.set_title(f"DES DR2 SExtractor theta_j2000 vs r-band image moments, {len(r)} galaxies "
                     f"(A > 1.5\", B/A < 0.5)\nverdict: {verdict}", fontsize=8)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        ncol = 6
        nrow = int(np.ceil(min(len(imgs), 24) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.0 * ncol, 2.1 * nrow), squeeze=False)
        for axx, (g, data, wcs), row in zip(axes.flat, imgs[:24], rows[:24]):
            ny, nx = data.shape
            yy, xx = np.mgrid[-0.5:ny, -0.5:nx]
            cra, cdec = wcs.wcs_pix2world(xx, yy, 0)
            ex, no = G.tangent_offsets(cra, cdec, g.ra, g.dec)
            axx.pcolormesh(ex, no, G.stretch(data, "asinh", 20, 99.7), cmap="gray_r", shading="flat", rasterized=True)
            for sgn, col, ls in ((1, "#2a6fdb", "-"), (-1, "#e07b39", "--")):
                # matplotlib angle is from +x (East here, x axis reversed below) toward +y (North): 90 - PA
                axx.add_patch(Ellipse((0, 0), 6 * g.a_image * PIX, 6 * g.b_image * PIX, angle=90 - sgn * g.theta_j2000,
                                      fill=False, ec=col, ls=ls, lw=1))
            h = a.size / 2
            axx.set_xlim(h, -h); axx.set_ylim(-h, h); axx.set_aspect("equal"); axx.set_xticks([]); axx.set_yticks([])
            axx.set_title(f"img {row['pa_image']:.0f}, cat {g.theta_j2000:.0f}", fontsize=6)
        for axx in axes.flat[len(imgs[:24]):]:
            axx.axis("off")
        fig.suptitle("r-band cutouts, North up, East left; 3x(A, B) ellipse with PA = +theta_j2000 (blue, solid) "
                     "and -theta_j2000 (orange, dashed)", fontsize=8)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


if __name__ == "__main__":
    main()
