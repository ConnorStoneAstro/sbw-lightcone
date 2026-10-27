"""The DC2 position-angle mirror: SNe were placed along +position_angle_true (East of North) while the
rendered galaxies, as the DP0.2 pipeline measures them, have PA = -position_angle_true. Left: SN offset
angles of elongated hosts relative to the pipeline moment PA, to its mirror image, and to the cosmoDC2 PA.
Right: coadd cutouts of large elongated hosts with their SN and both ellipses."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from dp02_fig import load
from dp02_gallery import (RateLimiter, _session, fetch_cutout, load_tokens, patch_table, read_image, stretch,
                          tangent_offsets)
from dp02_lib import COLORS, GAUSS_HALF_LIGHT, PIX, host_offsets, moment_shape
from sbw_host.core import sky_to_plane

Q_MAX = 0.5        # elongated: truth axis ratio of the dominant component and moment axis ratio below this
OFF_MIN = 1.0      # SN beyond this many R_e of the dominant component
N_EXAMPLES = 6
SHAPE_COLS = ["galaxy_id", "size_bulge_true", "size_minor_bulge_true", "size_disk_true", "size_minor_disk_true",
              "position_angle_true", "bulge_to_total_ratio_i"]


def fold(d):
    """Angle difference folded to [-90, 90) degrees (an axis, not a direction)."""
    return (np.asarray(d) + 90.0) % 180.0 - 90.0


def build(a):
    df = load(a.associations)
    s = pq.read_table(a.sample, columns=["sn_id", "host_galaxy", "host_ra", "host_dec"]).to_pandas(ignore_metadata=True)
    df = df.drop(columns=[c for c in ("host_ra", "host_dec") if c in df]).merge(s.drop(columns="host_galaxy"), on="sn_id")
    sh = pq.read_table(a.shapes, columns=SHAPE_COLS).to_pandas(ignore_metadata=True).rename(columns={"galaxy_id": "host_galaxy"})
    df = df.merge(sh, on="host_galaxy", how="inner")
    b = a.band
    mcols = [f"{b}_{k}" for k in ("ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF")]
    t = pq.read_table(a.candidates, columns=["objectId", "coord_ra", "coord_dec"] + mcols)
    t = t.filter(pc.is_in(t["objectId"], value_set=pa.array(df.host_objectId.unique(), pa.int64())))
    c = t.to_pandas(ignore_metadata=True).drop_duplicates("objectId")
    q, phi, l1, l2, valid = moment_shape(*(c[k].to_numpy(dtype=float) for k in mcols))
    # phi is from +x (West) toward +y (North); PA East of North of that axis
    c["q_mom"], c["pa_mom"] = q, np.degrees(np.arctan2(-np.cos(np.radians(phi)), np.sin(np.radians(phi))))
    c["a_mom"], c["b_mom"], c["mom_valid"] = np.sqrt(l1) * PIX, np.sqrt(l2) * PIX, valid
    df = df.merge(c[["objectId", "q_mom", "pa_mom", "a_mom", "b_mom", "mom_valid"]].rename(columns={"objectId": "host_objectId"}),
                  on="host_objectId", how="inner")
    off = host_offsets(df)
    bd = off["bulge_dominated"]
    comp = np.where(bd, "bulge", "disk")
    df["re"] = np.where(bd, df.size_bulge_true, df.size_disk_true)
    df["q_true"] = np.where(bd, df.size_minor_bulge_true / df.size_bulge_true, df.size_minor_disk_true / df.size_disk_true)
    df["comp"] = comp
    df["offset_re"] = off["own"]
    x, y = sky_to_plane(df.sn_ra.to_numpy(), df.sn_dec.to_numpy(), df.host_ra.to_numpy(), df.host_dec.to_numpy())
    df["psi"] = np.degrees(np.arctan2(-x, y))      # SN offset angle East of North
    df["pa_true"] = df.position_angle_true
    return df


def summary(sel):
    out = {"n": int(len(sel))}
    for name, ref in (("vs_moment_pa", sel.pa_mom), ("vs_mirrored_moment_pa", -sel.pa_mom), ("vs_cosmodc2_pa", sel.pa_true)):
        d = fold(sel.psi - ref)
        out[name] = {"mean_cos2": float(np.mean(np.cos(np.radians(2 * d)))), "median_abs_deg": float(np.median(np.abs(d)))}
    out["moment_pa_plus_cosmodc2_pa_median_abs_deg"] = float(np.median(np.abs(fold(sel.pa_mom + sel.pa_true))))
    out["moment_pa_minus_cosmodc2_pa_median_abs_deg"] = float(np.median(np.abs(fold(sel.pa_mom - sel.pa_true))))
    return out


def choose(df):
    """Large elongated hosts whose cosmoDC2 PA and measured PA differ by > 60 deg, SN within 15 deg of +PA
    at 1-3 R_e; the largest N_EXAMPLES, one per host."""
    m = ((df.q_true < 0.4) & (df.q_mom < 0.4) & df.mom_valid & (df.re > 1.5)
         & (np.abs(fold(df.pa_true - df.pa_mom)) > 60) & (np.abs(fold(df.psi - df.pa_true)) < 15)
         & (df.offset_re > 1.0) & (df.offset_re < 3.0))
    return df[m].sort_values("re", ascending=False).drop_duplicates("host_galaxy").head(N_EXAMPLES)


def ellipse(ax, e0, n0, a, b, pa, **kw):
    s = np.linspace(0, 2 * np.pi, 200)
    u, v = a * np.cos(s), b * np.sin(s)
    th = np.radians(pa)
    ax.plot(e0 + u * np.sin(th) + v * np.cos(th), n0 + u * np.cos(th) - v * np.sin(th), **kw)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--shapes", required=True)
    ap.add_argument("--band", required=True, choices=list("ugrizy"))
    ap.add_argument("--token-file", required=True)
    ap.add_argument("--cache-dir", default="cutout_cache")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = build(a)
    sel = df[(df.q_true < Q_MAX) & (df.q_mom < Q_MAX) & df.mom_valid & (df.offset_re > OFF_MIN)]
    summ = summary(sel)
    ex = choose(df)
    print(json.dumps({"elongated": summ, "examples": ex[["sn_id", "host_objectId", "re", "q_true", "q_mom", "pa_true",
                                                         "pa_mom", "psi", "offset_re"]].round(3).to_dict(orient="records")}),
          file=sys.stderr)

    fig = plt.figure(figsize=(12.5, 5.6))
    gs = fig.add_gridspec(2, 5, width_ratios=[1.5, 0.15, 1, 1, 1])
    ah = fig.add_subplot(gs[:, 0])
    edges = np.linspace(-90, 90, 37)
    for name, ref, color, ls in (("relative to the pipeline moment PA", sel.pa_mom, COLORS["dlr"], "-"),
                                 ("relative to the mirrored moment PA (-PA)", -sel.pa_mom, COLORS["sbw"], "-"),
                                 ("relative to the cosmoDC2 PA (+PA)", sel.pa_true, "black", ":")):
        d = fold(sel.psi - ref)
        ah.hist(d, bins=edges, histtype="step", color=color, ls=ls, lw=1.6, density=True,
                label=f"{name}\n<cos 2Δ> = {np.mean(np.cos(np.radians(2 * d))):.2f}")
    ah.axhline(1 / 180, color=COLORS["ref"], lw=0.8, ls="--")
    ah.set_xlabel("SN offset angle minus host PA [deg]"); ah.set_ylabel("density")
    ah.set_xlim(-90, 90); ah.set_xticks(np.arange(-90, 91, 30))
    ah.legend(fontsize=7, loc="upper left")
    ah.set_title(f"elongated hosts (q < {Q_MAX:g} in truth and moments), SN beyond {OFF_MIN:g} $R_e$: "
                 f"N = {len(sel):,}\nmedian |PA$_{{mom}}$ + PA$_{{cosmoDC2}}$| = "
                 f"{summ['moment_pa_plus_cosmodc2_pa_median_abs_deg']:.1f} deg", fontsize=8)

    tap_token, image_token = load_tokens(a.token_file)
    cache = Path(a.cache_dir)
    patches = patch_table(a.band, tap_token, cache)
    sess, lim = _session(image_token), RateLimiter(0.5)
    for k, (_, r) in enumerate(ex.iterrows()):
        ax = fig.add_subplot(gs[k // 3, 2 + k % 3])
        size = float(np.clip(7 * r.re, 15, 45))
        try:
            data, wcs = read_image(fetch_cutout(r.host_ra, r.host_dec, a.band, size, patches, sess, cache, lim))
        except Exception as e:   # noqa: BLE001 - leave the panel empty with the reason
            ax.text(0, 0, f"cutout failed:\n{type(e).__name__}", ha="center", fontsize=6); data = None
        if data is not None:
            ny, nx = data.shape
            yy, xx = np.mgrid[-0.5:ny, -0.5:nx]
            cra, cdec = wcs.wcs_pix2world(xx, yy, 0)
            e_, n_ = tangent_offsets(cra, cdec, r.host_ra, r.host_dec)
            ax.pcolormesh(e_, n_, stretch(data, "asinh", 20, 99.7), cmap="gray_r", shading="flat", rasterized=True)
        ellipse(ax, 0, 0, 2 * r.re, 2 * r.re * r.q_true, r.pa_true, color="black", ls=":", lw=1.4)
        ellipse(ax, 0, 0, 2 * GAUSS_HALF_LIGHT * r.a_mom, 2 * GAUSS_HALF_LIGHT * r.b_mom, r.pa_mom, color=COLORS["dlr"], lw=1.2)
        se, sn_ = tangent_offsets(r.sn_ra, r.sn_dec, r.host_ra, r.host_dec)
        ax.plot(se, sn_, "x", ms=9, mew=2, color="#d62728")
        h = size / 2
        ax.set_xlim(h, -h); ax.set_ylim(-h, h); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.grid(False)
        ax.set_title(f"SN {r.sn_id}: $R_e$ {r.re:.1f}\", q {r.q_true:.2f}\nPA$_{{cosmoDC2}}$ {r.pa_true:.0f}°, "
                     f"PA$_{{mom}}$ {r.pa_mom:.0f}°, P(true) {np.nan_to_num(r.p_sbw_true):.2f}", fontsize=7)
    fig.legend(handles=[plt.Line2D([], [], color="black", ls=":", label="cosmoDC2 ellipse (+PA), 2 $R_e$"),
                        plt.Line2D([], [], color=COLORS["dlr"], label="pipeline moment ellipse, 2 half-light"),
                        plt.Line2D([], [], ls="none", marker="x", color="#d62728", mew=2, label="SN")],
               loc="lower right", ncol=3, fontsize=7, frameon=False)
    fig.suptitle(f"DC2 placed SNe on the mirror image of the rendered galaxy ({a.band} band moments and deepCoadd cutouts, "
                 "north up, east left)", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(a.output); plt.close(fig)


if __name__ == "__main__":
    main()
