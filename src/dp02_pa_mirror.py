"""The DC2 position-angle mirror: SNe were placed along +position_angle_true (East of North) while the
rendered galaxies, as the DP0.2 pipeline measures them, have PA = -position_angle_true. Left: SN offset
angles of elongated hosts relative to the pipeline moment PA, to its mirror image, and to the cosmoDC2 PA.
Right: coadd cutouts of large elongated hosts with their SN and both ellipses. Page 2 (and ``--metrics``, a
json without cutouts): whether the SNe sit on the mirrored truth ellipse at the truth radii or on a
smaller / rounder one: SN radii along the major and minor truth axes, the rendered (measured) sizes against
the cosmoDC2 sizes, and the effective axis ratio of the SN distribution against the truth q."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
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
AXIS_DEG = 30.0    # SN "on the major axis" within this angle of +PA, "on the minor axis" within it of +PA + 90
QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9)
Q_BINS = np.array([0.05, 0.2, 0.3, 0.4, 0.5])
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
    t = pq.read_table(a.candidates, columns=["objectId", "coord_ra", "coord_dec", f"{b}_bdReB", f"{b}_bdReD"] + mcols)
    t = t.filter(pc.is_in(t["objectId"], value_set=pa.array(df.host_objectId.unique(), pa.int64())))
    c = t.to_pandas(ignore_metadata=True).drop_duplicates("objectId")
    q, phi, l1, l2, valid = moment_shape(*(c[k].to_numpy(dtype=float) for k in mcols))
    # phi is from +x (West) toward +y (North); PA East of North of that axis
    c["q_mom"], c["pa_mom"] = q, np.degrees(np.arctan2(-np.cos(np.radians(phi)), np.sin(np.radians(phi))))
    c["a_mom"], c["b_mom"], c["mom_valid"] = np.sqrt(l1) * PIX, np.sqrt(l2) * PIX, valid
    c["re_bulge_meas"], c["re_disk_meas"] = c[f"{b}_bdReB"] * PIX, c[f"{b}_bdReD"] * PIX   # determinant radii sqrt(ab)
    df = df.merge(c[["objectId", "q_mom", "pa_mom", "a_mom", "b_mom", "mom_valid", "re_bulge_meas", "re_disk_meas"]]
                  .rename(columns={"objectId": "host_objectId"}),
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
    df["sep"] = np.hypot(x, y)
    df["r_re"] = df.sep / df.re           # circular radius in units of the semi-major R_e
    d = np.radians(fold(df.psi - df.pa_true))
    df["u_re"], df["v_re"] = df.r_re * np.cos(d), df.r_re * np.sin(d)   # along / across the +PA (truth) major axis
    df["axis"] = np.where(np.abs(np.degrees(d)) < AXIS_DEG, "major",
                          np.where(np.abs(np.degrees(d)) > 90 - AXIS_DEG, "minor", "between"))
    return df


def summary(sel):
    out = {"n": int(len(sel))}
    for name, ref in (("vs_moment_pa", sel.pa_mom), ("vs_mirrored_moment_pa", -sel.pa_mom), ("vs_cosmodc2_pa", sel.pa_true)):
        d = fold(sel.psi - ref)
        out[name] = {"mean_cos2": float(np.mean(np.cos(np.radians(2 * d)))), "median_abs_deg": float(np.median(np.abs(d)))}
    out["moment_pa_plus_cosmodc2_pa_median_abs_deg"] = float(np.median(np.abs(fold(sel.pa_mom + sel.pa_true))))
    out["moment_pa_minus_cosmodc2_pa_median_abs_deg"] = float(np.median(np.abs(fold(sel.pa_mom - sel.pa_true))))
    return out


def quant(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return {f"q{int(100 * q)}": float(np.quantile(x, q)) for q in QUANTILES} | {"n": int(len(x))}


def radii(df):
    """SN radii of elongated hosts (truth q < Q_MAX, no offset cap): (a) circular and elliptical radii in R_e on
    and across the +PA axis, (b) measured over truth sizes, (c) the SN distribution's axis ratio against q."""
    e = df[df.q_true < Q_MAX]
    out = {"n_elongated": int(len(e)), "selection": f"dominant-component truth q < {Q_MAX:g}, all offsets",
           "axis_definition": f"major: SN within {AXIS_DEG:g} deg of +PA (the mirrored axis of the rendered galaxy); "
                              f"minor: within {AXIS_DEG:g} deg of +PA + 90"}
    def by_axis(e):
        ax = {}
        for k in ("major", "minor"):
            g = e[e.axis == k]
            ax[k] = {"n": int(len(g)), "r_over_re": quant(g.r_re), "elliptical_offset_re": quant(g.offset_re),
                     "median_q_true": float(g.q_true.median())}
        ax["minor_over_major_median_r"] = ax["minor"]["r_over_re"]["q50"] / ax["major"]["r_over_re"]["q50"]
        ax["minor_over_major_median_elliptical"] = (ax["minor"]["elliptical_offset_re"]["q50"]
                                                    / ax["major"]["elliptical_offset_re"]["q50"])
        return ax
    cap = float(np.ceil(df.sep.max()))   # the sample's separation cut (3 arcsec for clean3, 10 for clean10)
    out["a_radii_by_axis"] = by_axis(e)
    out["a_radii_by_axis_cap_safe"] = {"selection": f"6 R_e < {cap:g} arcsec, the sample's largest separation, "
                                                    "so the separation cut does not truncate the radii",
                                       **by_axis(e[6 * e.re < cap])}
    out["a_expectation"] = ("SB-weighted placement on the truth ellipse: equal elliptical-radius distributions on both "
                            "axes and minor/major circular radius = median q_true; a rounder placement ellipse gives "
                            "minor/major circular radius > q_true and elliptical radii larger on the minor axis")
    m = e[e.mom_valid & (e.q_mom < Q_MAX)]
    re_min = np.where(m.comp == "bulge", m.size_minor_bulge_true, m.size_minor_disk_true)
    meas_comp = np.where(m.comp == "bulge", m.re_bulge_meas, m.re_disk_meas)
    out["b_measured_over_truth_size"] = {
        "n": int(len(m)), "selection": "elongated in truth and in the moments (q < 0.5), valid moments",
        "moment_major_half_light_over_re_major": quant(GAUSS_HALF_LIGHT * m.a_mom / m.re),
        "moment_minor_half_light_over_re_minor": quant(GAUSS_HALF_LIGHT * m.b_mom / re_min),
        "q_moment_over_q_true": quant(m.q_mom / m.q_true),
        "cmodel_dominant_component_over_truth_sqrt_ab": quant(meas_comp / np.sqrt(m.re * re_min)),
        "by_component_and_size": {f"{c}_re_{lo:g}_{hi:g}": {
            "n": int(sel.sum()),
            "moment_major_over_re_major_q50": float(np.median(GAUSS_HALF_LIGHT * m.a_mom[sel] / m.re[sel])),
            "cmodel_over_truth_sqrt_ab_q50": float(np.nanmedian(meas_comp[sel] / np.sqrt(m.re * re_min)[sel])),
            "q_moment_over_q_true_q50": float(np.median(m.q_mom[sel] / m.q_true[sel]))}
            for c in ("bulge", "disk") for lo, hi in ((0, 0.5), (0.5, 1), (1, 99))
            for sel in [((m.comp == c) & (m.re >= lo) & (m.re < hi)).to_numpy()] if sel.sum() >= 20},
        "note": "moments: PSF-subtracted rms semi-axes x 1.1774 (Gaussian half-light; for a pure exponential this "
                "is 1.21 x the true half-light radius, more for n = 4); CModel: bdReB for "
                "bulge-dominated, bdReD for disk-dominated hosts, both determinant radii sqrt(ab), against the "
                "cosmoDC2 sqrt(a b) of the dominant component"}
    rows = []
    for lo, hi in zip(Q_BINS[:-1], Q_BINS[1:]):
        g = e[(e.q_true >= lo) & (e.q_true < hi) & np.isfinite(e.u_re)]
        if len(g) < 20:
            continue
        qm = float(g.q_true.median())
        q_mad = float(np.median(np.abs(g.v_re)) / np.median(np.abs(g.u_re)))
        g1 = g[g.r_re < 3]   # rms ratio without the far n = 4 tails
        q_rms = float(np.sqrt((g1.v_re ** 2).mean() / (g1.u_re ** 2).mean()))
        rows.append({"q_lo": float(lo), "q_hi": float(hi), "n": int(len(g)), "median_q_true": qm,
                     "q_sn_median_abs": q_mad, "q_sn_rms_r_lt_3re": q_rms,
                     "alpha_median_abs": float(np.log(q_mad) / np.log(qm)), "q_true_pow_0p4": qm ** 0.4})
    out["c_sn_axis_ratio"] = {"bins": rows, "definition": "q_SN = median |v| / median |u| of the SN offsets along "
                              "(u) and across (v) the +PA axis, stacked in R_e units; alpha from q_SN = q_true^alpha "
                              "(1 = SNe on the truth ellipse, 0 = round)"}
    return out


def choose(df):
    """Large elongated hosts whose cosmoDC2 PA and measured PA differ by > 60 deg, SN within 15 deg of +PA
    at 1.3-3.5 R_e; the largest N_EXAMPLES, one per host."""
    m = ((df.q_true < 0.4) & (df.q_mom < 0.4) & df.mom_valid & (df.re > 1.0)
         & (np.abs(fold(df.pa_true - df.pa_mom)) > 60) & (np.abs(fold(df.psi - df.pa_true)) < 15)
         & (df.offset_re > 1.3) & (df.offset_re < 3.5))
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
    ap.add_argument("--token-file", help="RSP token file (figure only)")
    ap.add_argument("--metrics", action="store_true", help="write the json of the selection and radius checks only")
    ap.add_argument("--cache-dir", default="cutout_cache")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = build(a)
    sel = df[(df.q_true < Q_MAX) & (df.q_mom < Q_MAX) & df.mom_valid & (df.offset_re > OFF_MIN)]
    summ = summary(sel)
    ex = choose(df)
    excols = ["sn_id", "host_objectId", "re", "q_true", "q_mom", "pa_true", "pa_mom", "psi", "offset_re", "r_re"]
    if a.metrics:
        pool = df[(df.q_true < 0.4) & (df.q_mom < 0.4) & df.mom_valid & (df.re > 1.0)
                  & (np.abs(fold(df.pa_true - df.pa_mom)) > 60) & (np.abs(fold(df.psi - df.pa_true)) < 15)]
        res = {"elongated_angles": summ, "radii": radii(df),
               "example_selection": {
                   "rule": "q < 0.4 in truth and moments, R_e > 1 arcsec, PAs differ > 60 deg, SN within 15 deg of "
                           "+PA, 1.3 < offset_re < 3.5, largest R_e first",
                   "examples": ex[excols].round(3).to_dict(orient="records"),
                   "offset_re_all_candidates_without_offset_cap": quant(pool.offset_re),
                   "frac_candidates_offset_re_lt_1p3": float((pool.offset_re < 1.3).mean()),
                   "frac_candidates_offset_re_gt_3p5": float((pool.offset_re > 3.5).mean())},
               "settings": {"band": a.band, "q_max": Q_MAX, "off_min": OFF_MIN}}
        with open(a.output, "w") as f:
            json.dump(res, f, indent=2, default=float)
        return
    if not a.token_file:
        raise SystemExit("--token-file is needed for the figure")
    print(json.dumps({"elongated": summ, "examples": ex[excols].round(3).to_dict(orient="records")}), file=sys.stderr)

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
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    with PdfPages(a.output) as pdf:
        pdf.savefig(fig); plt.close(fig)
        fig2 = radius_page(df)
        pdf.savefig(fig2); plt.close(fig2)


def radius_page(df):
    """Page 2: SN radii along the mirrored major and the minor truth axis, measured over truth sizes, and the
    effective axis ratio of the SN distribution."""
    r = radii(df)
    e = df[df.q_true < Q_MAX]
    fig, (aa, ab, ac, ad) = plt.subplots(1, 4, figsize=(15, 3.9))
    edges = np.linspace(0, 6, 49)
    for k, color in (("major", COLORS["sbw"]), ("minor", COLORS["dlr"])):
        g = e[e.axis == k]
        aa.hist(g.r_re, bins=edges, density=True, histtype="step", color=color, lw=1.5,
                label=f"{k} axis (N = {len(g):,}), median {g.r_re.median():.2f}")
        ab.hist(g.offset_re, bins=edges, density=True, histtype="step", color=color, lw=1.5,
                label=f"{k} axis, median {g.offset_re.median():.2f}")
    qm = float(e[e.axis == "minor"].q_true.median())
    aa.set_xlabel("SN radius / truth semi-major $R_e$ (circular)"); aa.set_ylabel("density")
    aa.legend(fontsize=7); aa.set_title(f"(a) SN radius, within {AXIS_DEG:g} deg of the +PA major or minor axis\n"
                                        f"truth ellipse would give minor/major = q ~ {qm:.2f}; "
                                        f"measured {r['a_radii_by_axis']['minor_over_major_median_r']:.2f}", fontsize=8)
    ab.set_xlabel("SN elliptical radius on the truth ellipse [$R_e$]"); ab.legend(fontsize=7)
    ab.set_title("(b) the same in elliptical radius (equal if SNe follow the truth ellipse)", fontsize=8)
    m = e[e.mom_valid & (e.q_mom < Q_MAX)]
    re_min = np.where(m.comp == "bulge", m.size_minor_bulge_true, m.size_minor_disk_true)
    le = np.linspace(-1.5, 1.5, 61)
    for y, lab, color in ((GAUSS_HALF_LIGHT * m.a_mom / m.re, "moment major / truth major", COLORS["sbw"]),
                          (GAUSS_HALF_LIGHT * m.b_mom / re_min, "moment minor / truth minor", COLORS["dlr"]),
                          (np.where(m.comp == "bulge", m.re_bulge_meas, m.re_disk_meas) / np.sqrt(m.re * re_min),
                           "CModel sqrt(ab) / truth sqrt(ab)", "black")):
        y = np.log10(np.asarray(y, float)); y = y[np.isfinite(y)]
        ac.hist(y, bins=le, density=True, histtype="step", color=color, lw=1.4, label=f"{lab}: median {10 ** np.median(y):.2f}")
    ac.axvline(0, color=COLORS["ref"], ls=":", lw=1)
    ac.set_xlabel("log10 measured / cosmoDC2 half-light size"); ac.legend(fontsize=7)
    ac.set_title(f"(c) rendered vs truth sizes (N = {len(m):,}, q < {Q_MAX:g} in both)", fontsize=8)
    b = pd.DataFrame(r["c_sn_axis_ratio"]["bins"])
    qq = np.linspace(0.05, 0.5, 50)
    ad.plot(qq, qq, color=COLORS["ref"], ls="--", lw=1, label="q_SN = q (truth ellipse)")
    ad.plot(qq, qq ** 0.4, color=COLORS["ref"], ls=":", lw=1, label="q_SN = q$^{0.4}$")
    if len(b):
        ad.plot(b.median_q_true, b.q_sn_median_abs, "o-", color=COLORS["sbw"], label="SNe, median |v| / median |u|")
        ad.plot(b.median_q_true, b.q_sn_rms_r_lt_3re, "s--", color=COLORS["sbw"], ms=4, lw=1, label="SNe, rms ratio (r < 3 $R_e$)")
    ad.set_xlabel("truth axis ratio q (dominant component)"); ad.set_ylabel("axis ratio of the SN distribution")
    ad.set_xlim(0, 0.55); ad.set_ylim(0, 1); ad.legend(fontsize=7)
    ad.set_title("(d) shape of the SN distribution about +PA", fontsize=8)
    fig.suptitle(f"DC2 SN radii on the mirrored (+PA) truth ellipse: {len(e):,} hosts with truth q < {Q_MAX:g}, all offsets",
                 fontsize=9)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    main()
