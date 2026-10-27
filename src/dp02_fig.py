"""Figures and the recovery metric from one DP0.2 association table."""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from dp02_lib import COLORS, DLR_MAX, HOSTLESS_CUT, calibration_curve, dlr_hostless, sbw_hostless, wilson

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5,
                     "legend.frameon": False, "figure.dpi": 120})


def load(path):
    t = pq.read_table(path)
    df = t.to_pandas(ignore_metadata=True)
    df.attrs = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if not k.startswith(b"pandas")}
    return df


def subtitle(df):
    a = df.attrs
    return f"{a.get('profile', '?')} | band {a.get('band', '?')} | {a.get('search_radius', '?')} | d_DLR<={a.get('dlr_max', '?')[1:]} ({a.get('dlr_radius', '?')}) | N={len(df)}"


def recovery(df, out, hostless_cut):
    n = len(df)
    cut = HOSTLESS_CUT[hostless_cut]
    hl = sbw_hostless(df, cut)
    kh = int((df.sbw_correct.to_numpy() & ~hl).sum())
    ph, lh, hh = wilson(kh, n)
    ks, kd = int(df.sbw_correct.sum()), int(df.dlr_correct.sum())
    ps, ls, hs = wilson(ks, n); pd_, ld, hd = wilson(kd, n)
    _, ece = calibration_curve(df)
    res = {"n": n, "sbw_recovery": ps, "sbw_recovery_ci68": [ls, hs],
           "dlr_recovery": pd_, "dlr_recovery_ci68": [ld, hd],
           "agreement": float(df.agree.mean()), "dlr_hostless_fraction": float(df.dlr_hostless.mean()),
           "sbw_no_host_fraction": float((df.sbw_objectId < 0).mean()),
           "true_host_in_candidates": float(df.true_in_cand.mean()),
           "median_p_sbw_true": float(np.nanmedian(df.p_sbw_true)),
           "expected_calibration_error": ece,
           "disagree_sbw_right": float((~df.agree & df.sbw_correct).mean()),
           "disagree_dlr_right": float((~df.agree & df.dlr_correct).mean()),
           "sbw_hostless_cut": cut, "sbw_hostless_fraction": float(hl.mean()),
           "sbw_recovery_with_hostless_cut": ph, "sbw_recovery_with_hostless_cut_ci68": [lh, hh],
           "settings": df.attrs}
    with open(out, "w") as f:
        json.dump(res, f, indent=2, default=float)


def plot_calibration(ax, df, label=None, color=COLORS["sbw"]):
    cc, ece = calibration_curve(df)
    cc = cc[cc.n > 0]
    ax.plot([0, 1], [0, 1], color=COLORS["ref"], lw=1, ls="--", zorder=1)
    ax.errorbar(cc.p_mean, cc.frac_correct, yerr=[cc.frac_correct - cc.lo, cc.hi - cc.frac_correct],
                fmt="o", ms=4, color=color, capsize=2, lw=1, label=label)
    ax.set_xlabel("top-host $P_{\\mathrm{SBW}}$"); ax.set_ylabel("fraction where top host is true host")
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02)
    ax.text(0.03, 0.95, f"ECE = {ece:.3f}", transform=ax.transAxes, va="top", fontsize=8)
    return cc, ece


def calibration(df, out):
    fig, (ax, ah) = plt.subplots(2, 1, figsize=(4.2, 5.2), height_ratios=[3, 1], sharex=True)
    plot_calibration(ax, df)
    ax.set_xlabel("")
    cc, _ = calibration_curve(df)
    ah.bar(cc.p_mean, cc.n, width=0.012, color=COLORS["sbw"], alpha=0.6)
    ah.set_yscale("log"); ah.set_ylabel("SNe per bin"); ah.set_xlabel("top-host $P_{\\mathrm{SBW}}$")
    fig.suptitle("SBW calibration on DC2 truth\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def agreement_matrix(df, out, hostless_cut):
    cut = HOSTLESS_CUT[hostless_cut]
    hl = sbw_hostless(df, cut)
    rows = [("SBW host = true", df.sbw_correct & ~hl), ("SBW host wrong", ~df.sbw_correct & ~hl),
            (f"SBW: no host\n(flb < {cut:g})", hl)]
    cols = [("DLR host = true", df.dlr_correct), ("DLR host wrong", ~df.dlr_correct & ~df.dlr_hostless),
            ("DLR: no host", df.dlr_hostless)]
    n = len(df)
    M = np.array([[int((r & c).sum()) for _, c in cols] for _, r in rows])
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    ax.imshow(M / n, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i, (rl, rm) in enumerate(rows):
        for j, (cl, cm) in enumerate(cols):
            sel = rm & cm
            txt = f"{M[i, j]:,}\n{100 * M[i, j] / n:.2f}%"
            if i == 1 and j == 1:
                same = int((sel & df.agree).sum())
                txt += f"\nsame wrong host: {same:,}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if M[i, j] / n > 0.5 else "black")
    ax.set_xticks(range(3), [c for c, _ in cols]); ax.set_yticks(range(3), [r for r, _ in rows])
    ax.grid(False)
    ax.set_title(f"SBW vs DLR vs truth (agree on host: {100 * df.agree.mean():.2f}%)\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def binned_fraction(x, flag, edges):
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
    ok = np.isfinite(x) & (x >= edges[0]) & (x < edges[-1])
    n = np.bincount(idx[ok], minlength=len(edges) - 1)
    k = np.bincount(idx[ok], weights=flag[ok].astype(float), minlength=len(edges) - 1)
    p, lo, hi = wilson(k, n)
    return 0.5 * (edges[1:] + edges[:-1]), p, lo, hi, n


def agreement_vs_ddlr(df, out):
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 4.8), height_ratios=[3, 1], sharex="col")
    for j, (col, lab) in enumerate([("d_dlr_sbw", "$d_{DLR}$ of the SBW host"), ("d_dlr_min", "smallest $d_{DLR}$ (DLR host)")]):
        x = df[col].to_numpy(dtype=float)
        edges = np.linspace(0, 6, 25)
        c, p, lo, hi, n = binned_fraction(x, df.agree.to_numpy(), edges)
        ax = axes[0, j]
        ax.fill_between(c, lo, hi, color=COLORS["agree"], alpha=0.2, step=None)
        ax.plot(c, p, color=COLORS["agree"], label="SBW host = DLR host")
        c2, p2, lo2, hi2, _ = binned_fraction(x, df.sbw_correct.to_numpy(), edges)
        ax.plot(c2, p2, color=COLORS["sbw"], label="SBW host = true host")
        c3, p3, lo3, hi3, _ = binned_fraction(x, df.dlr_correct.to_numpy(), edges)
        ax.plot(c3, p3, color=COLORS["dlr"], label="DLR host = true host")
        ax.set_ylim(0, 1.02); ax.set_ylabel("fraction")
        ax.axvline(float(df.attrs.get("dlr_max", "d4")[1:]), color=COLORS["ref"], ls=":", lw=1)
        if j == 0:
            ax.legend(fontsize=7, loc="lower left")
        axes[1, j].bar(c, n, width=edges[1] - edges[0], color=COLORS["ref"], alpha=0.5)
        axes[1, j].set_yscale("log"); axes[1, j].set_xlabel(lab); axes[1, j].set_ylabel("SNe")
    fig.suptitle("SBW-DLR agreement versus $d_{DLR}$\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


XSPEC = {
    "ddlr": ("d_dlr_true", "$d_{DLR}$ of the true host", np.linspace(0, 6, 25), False),
    "separation": ("sep_true", "SN-host separation [arcsec]", np.geomspace(0.05, 30, 25), True),
    "host_mag": ("host_mag_r", "true host $r$ [mag]", np.linspace(16, 29, 27), False),
    "redshift": ("sn_redshift", "SN redshift", np.linspace(0, 1.5, 31), False),
}


def recovery_vs(df, out, xkey):
    col, lab, edges, logx = XSPEC[xkey]
    x = df[col].to_numpy(dtype=float)
    fig, (ax, ah) = plt.subplots(2, 1, figsize=(4.8, 4.6), height_ratios=[3, 1], sharex=True)
    for flag, name, color in [(df.sbw_correct, "SBW", COLORS["sbw"]), (df.dlr_correct, "DLR", COLORS["dlr"])]:
        c, p, lo, hi, n = binned_fraction(x, flag.to_numpy(), edges)
        ax.fill_between(c, lo, hi, color=color, alpha=0.15)
        ax.plot(c, p, color=color, label=f"{name} host = true host")
    c, p, lo, hi, n = binned_fraction(x, df.agree.to_numpy(), edges)
    ax.plot(c, p, color=COLORS["agree"], ls="--", lw=1, label="SBW host = DLR host")
    ax.set_ylim(0, 1.02); ax.set_ylabel("fraction"); ax.legend(fontsize=7, loc="lower left")
    ah.bar(c, n, width=np.diff(edges), color=COLORS["ref"], alpha=0.5)
    ah.set_yscale("log"); ah.set_ylabel("SNe"); ah.set_xlabel(lab)
    if logx:
        ax.set_xscale("log")
    fig.suptitle("Host recovery on DC2 truth\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def p_true_host(df, out):
    p = df.p_sbw_true.to_numpy(dtype=float)
    ok = np.isfinite(p) & (p > 0)
    lp = np.log10(np.clip(p, 1e-12, 1))
    edges = np.linspace(-12, 0, 49)
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.hist(lp[ok & df.sbw_correct.to_numpy()], bins=edges, color=COLORS["sbw"], alpha=0.7, label="true host ranked first")
    ax.hist(lp[ok & ~df.sbw_correct.to_numpy()], bins=edges, color=COLORS["dlr"], alpha=0.7, label="true host not first")
    ax.set_yscale("log"); ax.set_xlabel("$\\log_{10} P_{\\mathrm{SBW}}$(true host)"); ax.set_ylabel("SNe")
    nz = int((np.isfinite(p) & (p <= 0)).sum()); nmiss = int((~np.isfinite(p)).sum())
    ax.text(0.02, 0.95, f"P = 0: {nz:,}\ntrue host not a candidate: {nmiss:,}", transform=ax.transAxes, va="top", fontsize=7)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.0, 0.8))
    fig.suptitle("Probability SBW gives the true host\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


# ---------------------------------------------------------------- hostless indicator
FLB_GRID = np.geomspace(1e-5, 1, 121)
DLR_GRID = np.linspace(0, 12, 121)


def rejected(flag):
    flag = np.asarray(flag)
    return wilson(int(flag.sum()), len(flag))


def sweep(df, rnd, rule, grid):
    """Real-kept and random-rejected fractions (with Wilson bands) for every threshold in ``grid``."""
    real = np.array([rejected(~rule(df, t)) for t in grid])
    rand = np.array([rejected(rule(rnd, t)) for t in grid])
    return real, rand


RULES = {"flb_best": lambda d, t: sbw_hostless(d, t, "flb_best"),
         "flb_max": lambda d, t: sbw_hostless(d, t, "flb_max"),
         "d_dlr_min": dlr_hostless}


def hostless_metrics(df, rnd, out):
    def point(rule, t):
        kept = ~RULES[rule](df, t); rej = RULES[rule](rnd, t)
        pk, lk, hk = wilson(int(kept.sum()), len(df)); pr, lr, hr = wilson(int(rej.sum()), len(rnd))
        on_gal = (~rej) & (rnd.p_sbw_max.to_numpy(dtype=float) > 0.5)
        return {"threshold": t, "real_kept": pk, "real_kept_ci68": [lk, hk],
                "real_kept_and_sbw_correct": float((kept & df.sbw_correct.to_numpy()).mean()),
                "random_rejected": pr, "random_rejected_ci68": [lr, hr],
                "random_kept_with_p_sbw_gt_0p5": float(on_gal.mean())}
    p = rnd.p_sbw_max.to_numpy(dtype=float)
    res = {"n_real": len(df), "n_random": len(rnd),
           "random_with_no_candidate": float((rnd.n_cand_used == 0).mean()),
           "random_with_no_sbw_host": float((rnd.sbw_objectId < 0).mean()),
           "random_p_sbw_gt_0p5": float((p > 0.5).mean()),
           "random_p_sbw_gt_0p9": float((p > 0.9).mean()),
           "real_with_no_sbw_host": float((df.sbw_objectId < 0).mean()),
           "sbw_flb_best": [point("flb_best", t) for t in (0.001, 0.01, 0.05, 0.1)],
           "sbw_flb_max": [point("flb_max", t) for t in (0.001, 0.01, 0.05, 0.1)],
           "dlr_d_dlr_min": [point("d_dlr_min", t) for t in (3.0, 4.0, 5.0)],
           "rule": "real kept = not hostless; SBW hostless = no SBW host or flb < t; DLR hostless = no candidate with d_DLR <= t",
           "settings": df.attrs}
    with open(out, "w") as f:
        json.dump(res, f, indent=2, default=float)


def band(ax, x, w, color, ls="-", label=None):
    ax.fill_between(x, w[:, 1], w[:, 2], color=color, alpha=0.15, lw=0)
    ax.plot(x, w[:, 0], color=color, ls=ls, label=label)


def hostless_indicator(df, rnd, out, hostless_cut, dlr_max):
    cut, dmax = HOSTLESS_CUT[hostless_cut], DLR_MAX[dlr_max]
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.0))
    (aa, ab), (ac, ad) = axes
    # (a) flb_best distributions
    edges = np.geomspace(1e-4, 1, 41)
    lf = lambda d: np.clip(d.flb_best.to_numpy(dtype=float), edges[0] * 1.0001, 1)
    real_ok = df.sbw_objectId.to_numpy() >= 0
    for sel, d, lab, col, kw in [(real_ok & df.sbw_correct.to_numpy(), df, "real SN, SBW host = true", COLORS["sbw"], {"alpha": 0.6}),
                                 (real_ok & ~df.sbw_correct.to_numpy(), df, "real SN, SBW host wrong", COLORS["dlr"], {"alpha": 0.6}),
                                 (rnd.sbw_objectId.to_numpy() >= 0, rnd, "random position", COLORS["agree"], {"histtype": "step", "lw": 1.5})]:
        aa.hist(lf(d)[sel], bins=edges, weights=np.full(sel.sum(), 1 / len(d)), color=col, label=lab, **kw)
    aa.set_xscale("log"); aa.set_yscale("log")
    aa.set_xlabel("frac_light_beyond of the SBW host (flb_best; < $10^{-4}$ in first bin)")
    aa.set_ylabel("fraction of population per bin")
    aa.axvline(cut, color=COLORS["ref"], ls=":", lw=1)
    aa.text(0.02, 0.97, f"no SBW host: real {100 * (~real_ok).mean():.2f}%, random {100 * (rnd.sbw_objectId < 0).mean():.2f}%",
            transform=aa.transAxes, va="top", fontsize=7)
    aa.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0, 0.92))
    aa.set_title("(a) SBW hostless statistic", fontsize=9, loc="left")
    # (b) SBW rejection vs threshold
    fb_real, fb_rand = sweep(df, rnd, RULES["flb_best"], FLB_GRID)
    band(ab, FLB_GRID, 1 - fb_real[:, [0, 2, 1]], COLORS["sbw"], label="real SNe rejected")
    band(ab, FLB_GRID, fb_rand, COLORS["agree"], ls="--", label="random positions rejected")
    ab.set_xscale("log"); ab.set_ylim(0, 1.02); ab.axvline(cut, color=COLORS["ref"], ls=":", lw=1)
    ab.set_xlabel("threshold t (hostless if flb_best < t or no SBW host)"); ab.set_ylabel("fraction rejected as hostless")
    ab.legend(fontsize=7, loc="upper left"); ab.set_title("(b) SBW: flb_best < t", fontsize=9, loc="left")
    # (c) DLR rejection vs threshold
    dd_real, dd_rand = sweep(df, rnd, RULES["d_dlr_min"], DLR_GRID)
    band(ac, DLR_GRID, 1 - dd_real[:, [0, 2, 1]], COLORS["dlr"], label="real SNe rejected")
    band(ac, DLR_GRID, dd_rand, COLORS["agree"], ls="--", label="random positions rejected")
    ac.set_ylim(0, 1.02); ac.axvline(dmax, color=COLORS["ref"], ls=":", lw=1)
    ac.set_xlabel("threshold t (hostless if $d_{DLR,\\min}$ > t)"); ac.set_ylabel("fraction rejected as hostless")
    ac.legend(fontsize=7, loc="center right"); ac.set_title("(c) DLR: $d_{DLR,\\min}$ > t", fontsize=9, loc="left")
    # (d) ROC
    fm_real, fm_rand = sweep(df, rnd, RULES["flb_max"], FLB_GRID)
    for (re_, ra_), lab, col, ls in [((fb_real, fb_rand), "SBW, flb_best", COLORS["sbw"], "-"),
                                     ((fm_real, fm_rand), "SBW, flb_max", COLORS["sbw"], "--"),
                                     ((dd_real, dd_rand), "DLR, $d_{DLR,\\min}$", COLORS["dlr"], "-")]:
        ad.plot(re_[:, 0], ra_[:, 0], color=col, ls=ls, label=lab)
    for rule, t, col, lab in [("flb_best", cut, COLORS["sbw"], f"flb_best = {cut:g}"),
                              ("d_dlr_min", dmax, COLORS["dlr"], f"$d_{{DLR}}$ = {dmax:g}")]:
        k = (~RULES[rule](df, t)).mean(); r = RULES[rule](rnd, t).mean()
        ad.plot(k, r, "o", color=col, ms=6, mec="black", mew=0.6, label=f"{lab}: kept {100 * k:.1f}%, rejected {100 * r:.1f}%")
    ad.set_xlabel("fraction of real SNe kept (have a host)"); ad.set_ylabel("fraction of random positions rejected")
    ad.set_xlim(0.8, 1.002); ad.set_ylim(0, 1.02)
    ad.legend(fontsize=7, loc="lower left"); ad.set_title("(d) hostless indicators compared", fontsize=9, loc="left")
    fig.suptitle(f"Hostless indicators on DC2: {len(df):,} real SNe vs {len(rnd):,} random positions\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind", choices=["recovery", "calibration", "agreement_matrix", "agreement_vs_ddlr", "recovery_vs", "p_true_host",
                                    "hostless_indicator", "hostless_metrics"])
    ap.add_argument("--associations", required=True)
    ap.add_argument("--x", choices=sorted(XSPEC))
    ap.add_argument("--random", help="association table of random positions (hostless kinds)")
    ap.add_argument("--hostless-cut", choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--dlr-max", choices=sorted(DLR_MAX))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = load(a.associations)
    if a.kind == "recovery":
        recovery(df, a.output, a.hostless_cut)
    elif a.kind == "calibration":
        calibration(df, a.output)
    elif a.kind == "agreement_matrix":
        agreement_matrix(df, a.output, a.hostless_cut)
    elif a.kind == "agreement_vs_ddlr":
        agreement_vs_ddlr(df, a.output)
    elif a.kind == "recovery_vs":
        recovery_vs(df, a.output, a.x)
    elif a.kind == "p_true_host":
        p_true_host(df, a.output)
    elif a.kind == "hostless_indicator":
        hostless_indicator(df, load(a.random), a.output, a.hostless_cut, a.dlr_max)
    elif a.kind == "hostless_metrics":
        hostless_metrics(df, load(a.random), a.output)


if __name__ == "__main__":
    main()
