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

from dp02_lib import (CAL_EDGES, COLORS, DLR_MAX, HOSTLESS_CUT, P_TRUE, calibration_curve, dlr_cols, dlr_hostless, matrix_cells,
                      p_true, prob_matrix_cells, prob_rows, random_cells, sbw_hostless, wilson)

TIER = {"dp02": "DC2", "ou24": "OU2024", "dp2": "DP2"}
LABEL = {"tier": "DC2"}   # set from --tier; the figure titles name the simulation

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
    return (f"{a.get('profile', '?')} | band {a.get('band', '?')} | {a.get('search_radius', '?')} | d_DLR<={a.get('dlr_max', '?')[1:]} "
            f"({a.get('dlr_radius', '?')}) | quality {a.get('candidate_quality', '?')} | N={len(df)}")


def recovery(df, out, hostless_cut, rnd, p_threshold):
    n = len(df)
    cut, t = HOSTLESS_CUT[hostless_cut], P_TRUE[p_threshold]
    pt = p_true(df)
    bands = {k: float(v.mean()) for k, v in prob_rows(df, cut, t).items()}
    hl = sbw_hostless(df, cut)
    kh = int((df.sbw_correct.to_numpy() & ~hl).sum())
    ph, lh, hh = wilson(kh, n)
    ks, kd = int(df.sbw_correct.sum()), int(df.dlr_correct.sum())
    ps, ls, hs = wilson(ks, n); pd_, ld, hd = wilson(kd, n)
    _, ece, n_cal = calibration_curve(df, cut)
    res = {"n": n,
           "mean_p_sbw_true": float(pt.mean()), "mean_p_sbw_true_se": float(pt.std(ddof=1) / np.sqrt(n)),
           "mean_p_sbw_true_with_hostless_cut": float(np.where(hl, 0.0, pt).mean()),
           "mean_p_sbw_true_note": "expected recovery if each SN's host were drawn from the SBW probabilities; "
                                   "P = 0 where the true host is not a surviving candidate",
           "p_true_threshold": t, "p_true_bands": bands,
           "sbw_argmax_recovery": ps, "sbw_argmax_recovery_ci68": [ls, hs],
           "dlr_recovery": pd_, "dlr_recovery_ci68": [ld, hd],
           "agreement": float(df.agree.mean()), "dlr_hostless_fraction": float(df.dlr_hostless.mean()),
           "sbw_no_host_fraction": float((df.sbw_objectId < 0).mean()),
           "true_host_in_candidates": float(df.true_in_cand.mean()),
           "median_p_sbw_true": float(np.nanmedian(df.p_sbw_true)),
           "expected_calibration_error": ece, "calibration_n_sne": n_cal,
           "expected_calibration_error_note": "candidate-level: sum over P bins of |true hosts - summed P| per SN, on SNes SBW keeps "
                                              "as hosted (flb_best >= cut) whose true host is a surviving candidate",
           "disagree_sbw_right": float((~df.agree & df.sbw_correct).mean()),
           "disagree_dlr_right": float((~df.agree & df.dlr_correct).mean()),
           "sbw_hostless_cut": cut, "sbw_hostless_fraction": float(hl.mean()),
           "sbw_recovery_with_hostless_cut": ph, "sbw_recovery_with_hostless_cut_ci68": [lh, hh],
           "agreement_matrix": prob_matrix_cells(df, cut, t),
           "argmax_agreement_matrix": matrix_cells(df, cut),
           "random_matrix": random_cells(rnd, cut),
           "settings": df.attrs}
    with open(out, "w") as f:
        json.dump(res, f, indent=2, default=float)


def plot_calibration(ax, df, cut, label=None, color=COLORS["sbw"], log=True):
    cc, ece, nsn = calibration_curve(df, cut)
    cc = cc[cc.n_cand > 0]
    lim = (CAL_EDGES[0] * 0.8, 1.05)
    ax.plot(lim, lim, color=COLORS["ref"], lw=1, ls="--", zorder=1)
    ax.errorbar(cc.p_mean, cc.frac_true, yerr=[np.clip(cc.frac_true - cc.lo, 0, None), np.clip(cc.hi - cc.frac_true, 0, None)],
                fmt="o", ms=4, color=color, capsize=2, lw=1, label=label)
    if log:
        ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel("$P_{\\mathrm{SBW}}$ of a candidate"); ax.set_ylabel("fraction of candidates that are the true host")
    ax.text(0.03, 0.95, f"calibration error {ece:.3f} per SN\nN = {nsn:,} hosted SNe", transform=ax.transAxes, va="top", fontsize=8)
    return cc, ece


def calibration(df, out, hostless_cut):
    cut = HOSTLESS_CUT[hostless_cut]
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 5.6), height_ratios=[3, 1], sharex="col")
    (ax, az), (ah, ahz) = axes
    cc, ece = plot_calibration(ax, df, cut)
    plot_calibration(az, df, cut, log=False)
    az.set_xlim(0.45, 1.02); az.set_ylim(0.45, 1.02); az.set_ylabel(""); az.set_title("P > 0.5, linear", fontsize=8)
    ax.set_xlabel(""); az.set_xlabel("")
    for axis in (ah, ahz):
        axis.bar(cc.p_mean, cc.n_cand, width=0.08 * cc.p_mean if axis is ah else 0.012, color=COLORS["ref"], alpha=0.5, label="candidates")
        axis.bar(cc.p_mean, cc.n_true, width=0.05 * cc.p_mean if axis is ah else 0.008, color=COLORS["sbw"], alpha=0.8, label="true hosts")
        axis.set_yscale("log"); axis.set_xlabel("$P_{\\mathrm{SBW}}$ of a candidate")
    ah.set_xscale("log"); ah.set_ylabel("per bin"); ah.legend(fontsize=7, loc="upper left")
    fig.suptitle(f"SBW calibration on {LABEL['tier']} truth, every candidate with P > {CAL_EDGES[0]:g}\n"
                 f"SNe SBW keeps as hosted (flb_best >= {cut:g}) whose true host is a surviving candidate\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def agreement_matrix(df, out, hostless_cut, rnd, p_threshold, survivors_only=True):
    """Real-SN panel restricted (survivors_only) to SNe whose true host passed the catalogue (cat1) and
    object-table (cat2) cuts for both methods, so every cell is an association outcome; the random-position
    panel is unchanged."""
    cut, t = HOSTLESS_CUT[hostless_cut], P_TRUE[p_threshold]
    full = len(df)
    if survivors_only:
        survived = (df["true_host_status"].to_numpy() == "kept") & np.isfinite(df["d_dlr_true"].to_numpy(dtype=float))
        df = df[survived]
    pr, dc = prob_rows(df, cut, t), dlr_cols(df)
    rows = [("$P_{\\mathrm{SBW}}$(true) > 0.5", pr["p_gt_0p5"]),
            (f"$P_{{\\mathrm{{SBW}}}}$(true)\nin [{t:g}, 0.5]", pr["p_t_to_0p5"]),
            (f"$P_{{\\mathrm{{SBW}}}}$(true) < {t:g}\n(excluded)", pr["p_lt_t"]),
            (f"SBW: no host\n(flb < {cut:g})", pr["sbw_hostless"])]
    cols = [("DLR host = true", dc["dlr_true"]), ("DLR host wrong", dc["dlr_wrong"]), ("DLR: no host", dc["dlr_hostless"])]
    n = len(df)
    M = np.array([[int((r & c).sum()) for _, c in cols] for _, r in rows])
    fig, (ax, ar) = plt.subplots(1, 2, figsize=(9.8, 5.0), width_ratios=[3, 2.2])
    ax.imshow(M / n, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(cols)):
            ax.text(j, i, f"{M[i, j]:,}\n{100 * M[i, j] / n:.2f}%", ha="center", va="center", fontsize=8,
                    color="white" if M[i, j] / n > 0.5 else "black")
    ax.set_xticks(range(3), [c for c, _ in cols]); ax.set_yticks(range(len(rows)), [r for r, _ in rows])
    ax.grid(False)
    pt = np.where(pr["sbw_hostless"], 0.0, p_true(df))
    ax.set_title(f"real SNe whose true host passed cat1 + cat2: N = {n:,} of {full:,} ({100 * n / full:.1f}%)\n"
                 f"SBW most probable host = DLR host {100 * df.agree.mean():.2f}%", fontsize=8)
    ax.set_xlabel(f"mean $P_{{\\mathrm{{SBW}}}}$(true) with the hostless cut {100 * pt.mean():.2f}%; "
                  f"DLR recovery {100 * df.dlr_correct.mean():.2f}%", fontsize=8)
    # random positions: no host exists, so hostless is the correct answer
    hr = sbw_hostless(rnd, cut); dr = rnd["dlr_objectId"].to_numpy() < 0
    rrows = [("SBW: no host\n(correct)", hr), ("SBW assigns host\n(wrong)", ~hr)]
    rcols = [("DLR: no host\n(correct)", dr), ("DLR assigns host\n(wrong)", ~dr)]
    nr = len(rnd)
    R = np.array([[int((r & c).sum()) for _, c in rcols] for _, r in rrows])
    ar.imshow(R / nr, cmap="Greys", vmin=0, vmax=1.6, aspect="auto")
    for i in range(2):
        for j in range(2):
            ar.text(j, i, f"{R[i, j]:,}\n{100 * R[i, j] / nr:.2f}%", ha="center", va="center", fontsize=8)
    ar.set_xticks(range(2), [c for c, _ in rcols]); ar.set_yticks(range(2), [r for r, _ in rrows])
    ar.grid(False)
    ar.set_title(f"random positions (N = {nr:,})", fontsize=8)
    ar.set_xlabel(f"correctly rejected: SBW {100 * hr.mean():.2f}%, DLR {100 * dr.mean():.2f}%", fontsize=8)
    fig.suptitle("SBW probability of the true host vs DLR on SNe whose true host survived the catalogue and object-table cuts; "
                 "same hostless cuts on random positions\n"
                 + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def binned_fraction(x, flag, edges):
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
    ok = np.isfinite(x) & (x >= edges[0]) & (x < edges[-1])
    n = np.bincount(idx[ok], minlength=len(edges) - 1)
    k = np.bincount(idx[ok], weights=flag[ok].astype(float), minlength=len(edges) - 1)
    p, lo, hi = wilson(k, n)
    return 0.5 * (edges[1:] + edges[:-1]), p, lo, hi, n


def binned_mean(x, y, edges):
    """Bin centres, mean of y and its standard error per bin."""
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
    ok = np.isfinite(x) & (x >= edges[0]) & (x < edges[-1])
    n = np.bincount(idx[ok], minlength=len(edges) - 1)
    s1 = np.bincount(idx[ok], weights=y[ok], minlength=len(edges) - 1)
    s2 = np.bincount(idx[ok], weights=y[ok] ** 2, minlength=len(edges) - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = s1 / n
        se = np.sqrt(np.clip(s2 / n - m ** 2, 0, None) / np.maximum(n - 1, 1))
    return 0.5 * (edges[1:] + edges[:-1]), m, se


def agreement_vs_ddlr(df, out):
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 4.8), height_ratios=[3, 1], sharex="col")
    for j, (col, lab) in enumerate([("d_dlr_sbw", "$d_{DLR}$ of the SBW most probable host"), ("d_dlr_min", "smallest $d_{DLR}$ (DLR host)")]):
        x = df[col].to_numpy(dtype=float)
        edges = np.linspace(0, 6, 25)
        c, p, lo, hi, n = binned_fraction(x, df.agree.to_numpy(), edges)
        ax = axes[0, j]
        ax.fill_between(c, lo, hi, color=COLORS["agree"], alpha=0.2, step=None)
        ax.plot(c, p, color=COLORS["agree"], label="SBW most probable host = DLR host")
        c2, p2, _ = binned_mean(x, p_true(df), edges)
        ax.plot(c2, p2, color=COLORS["sbw"], label="mean $P_{\\mathrm{SBW}}$(true host)")
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
    if xkey == "redshift" and np.nanmax(x) > edges[-1]:   # OU2024 transients reach z ~ 3
        edges = np.linspace(0, 3, 31)
    fig, (ax, ah) = plt.subplots(2, 1, figsize=(4.8, 4.6), height_ratios=[3, 1], sharex=True)
    c, m, se = binned_mean(x, p_true(df), edges)
    ax.fill_between(c, m - se, m + se, color=COLORS["sbw"], alpha=0.15)
    ax.plot(c, m, color=COLORS["sbw"], label="mean $P_{\\mathrm{SBW}}$(true host)")
    c, p, lo, hi, n = binned_fraction(x, df.sbw_correct.to_numpy(), edges)
    ax.plot(c, p, color=COLORS["sbw"], ls="--", lw=1, label="SBW most probable host = true host")
    c, p, lo, hi, n = binned_fraction(x, df.dlr_correct.to_numpy(), edges)
    ax.fill_between(c, lo, hi, color=COLORS["dlr"], alpha=0.15)
    ax.plot(c, p, color=COLORS["dlr"], label="DLR host = true host")
    c, p, lo, hi, n = binned_fraction(x, df.agree.to_numpy(), edges)
    ax.plot(c, p, color=COLORS["agree"], ls=":", lw=1, label="SBW most probable host = DLR host")
    ax.set_ylim(0, 1.02); ax.set_ylabel("fraction"); ax.legend(fontsize=7, loc="lower left")
    ah.bar(c, n, width=np.diff(edges), color=COLORS["ref"], alpha=0.5)
    ah.set_yscale("log"); ah.set_ylabel("SNe"); ah.set_xlabel(lab)
    if logx:
        ax.set_xscale("log")
    fig.suptitle(f"Host recovery on {LABEL['tier']} truth\n" + subtitle(df), fontsize=8)
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
    ax.text(0.02, 0.95, f"P = 0: {nz:,}\ntrue host not a surviving candidate: {nmiss:,}", transform=ax.transAxes, va="top", fontsize=7)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.0, 0.8))
    fig.suptitle("Probability SBW gives the true host\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


# ---------------------------------------------------------------- hostless indicator
FLB_GRID = np.geomspace(1e-5, 1, 121)
DLR_GRID = np.linspace(0, 12, 121)
# Thresholds frozen from the DP0.2 baseline (27 Oct 2026: Youden-J optima of dp02_hostless_metrics and the knee
# of dp02_credible_set_metrics), next to the defaults; every tier reports itself at these values, so OU2024
# is a held-out test of whether the DP0.2 optima transfer.
DP02_OPTIMA = {"flb_best": 0.0215, "d_dlr_min": 2.7, "credible_set_t": 0.034}
TRANSFER = {"flb_best": (0.01, DP02_OPTIMA["flb_best"], 0.05), "d_dlr_min": (DP02_OPTIMA["d_dlr_min"], 4.0)}


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


def youden(df, rnd, rule, grid):
    """Youden's J = real kept + random rejected - 1 for every threshold, and the J-optimal threshold."""
    real, rand = sweep(df, rnd, RULES[rule], grid)
    j = real[:, 0] + rand[:, 0] - 1
    k = int(np.argmax(j))
    return j, {"t_opt": float(grid[k]), "J_max": float(j[k]), "real_kept": float(real[k, 0]),
               "random_rejected": float(rand[k, 0])}


def transfer(df, rnd):
    """Real kept, random rejected and J at the frozen DP0.2 thresholds and the defaults, next to this table's own
    J optimum and the J lost by using the DP0.2 value instead."""
    out = {"dp02_optima": DP02_OPTIMA}
    for rule, ts in TRANSFER.items():
        _, own = youden(df, rnd, rule, FLB_GRID if rule != "d_dlr_min" else DLR_GRID)
        pts = []
        for t in ts:
            k, r = float((~RULES[rule](df, t)).mean()), float(RULES[rule](rnd, t).mean())
            pts.append({"threshold": t, "real_kept": k, "random_rejected": r, "J": k + r - 1,
                        "J_minus_own_max": k + r - 1 - own["J_max"]})
        out[rule] = {"own_optimum": own, "at": pts}
    return out


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
           "youden_J": {"definition": "J = real kept + random rejected - 1, maximised over the threshold grid",
                        "sbw_flb_best": youden(df, rnd, "flb_best", FLB_GRID)[1],
                        "sbw_flb_max": youden(df, rnd, "flb_max", FLB_GRID)[1],
                        "dlr_d_dlr_min": youden(df, rnd, "d_dlr_min", DLR_GRID)[1],
                        "flb_grid": [float(FLB_GRID[0]), float(FLB_GRID[-1]), len(FLB_GRID), "log"],
                        "d_dlr_grid": [float(DLR_GRID[0]), float(DLR_GRID[-1]), len(DLR_GRID), "linear"]},
           "dp02_optima_transfer": transfer(df, rnd),
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
    if (df.host_objectId.to_numpy() >= 0).any():
        real = [(real_ok & df.sbw_correct.to_numpy(), df, "real SN, most probable host = true", COLORS["sbw"], {"alpha": 0.6}),
                (real_ok & ~df.sbw_correct.to_numpy(), df, "real SN, most probable host wrong", COLORS["dlr"], {"alpha": 0.6})]
    else:   # no truth (DP2)
        real = [(real_ok, df, "real transient", COLORS["sbw"], {"alpha": 0.6})]
    for sel, d, lab, col, kw in real + [(rnd.sbw_objectId.to_numpy() >= 0, rnd, "random position", COLORS["agree"], {"histtype": "step", "lw": 1.5})]:
        aa.hist(lf(d)[sel], bins=edges, weights=np.full(sel.sum(), 1 / len(d)), color=col, label=lab, **kw)
    aa.set_xscale("log"); aa.set_yscale("log")
    aa.set_xlabel("frac_light_beyond of the SBW most probable host (flb_best; < $10^{-4}$ in first bin)")
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
    j, jo = youden(df, rnd, "flb_best", FLB_GRID)
    ab.plot(FLB_GRID, j, color="black", lw=1, label="Youden J = kept + rejected - 1")
    ab.plot(jo["t_opt"], jo["J_max"], "v", color="black", ms=6, label=f"max J = {jo['J_max']:.3f} at t = {jo['t_opt']:.3g}")
    ab.set_xscale("log"); ab.set_ylim(0, 1.02); ab.axvline(cut, color=COLORS["ref"], ls=":", lw=1)
    ab.set_xlabel("threshold t (hostless if flb_best < t or no SBW host)"); ab.set_ylabel("fraction rejected as hostless")
    ab.legend(fontsize=7, loc="upper left"); ab.set_title("(b) SBW: flb_best < t", fontsize=9, loc="left")
    # (c) DLR rejection vs threshold
    dd_real, dd_rand = sweep(df, rnd, RULES["d_dlr_min"], DLR_GRID)
    band(ac, DLR_GRID, 1 - dd_real[:, [0, 2, 1]], COLORS["dlr"], label="real SNe rejected")
    band(ac, DLR_GRID, dd_rand, COLORS["agree"], ls="--", label="random positions rejected")
    j, jo = youden(df, rnd, "d_dlr_min", DLR_GRID)
    ac.plot(DLR_GRID, j, color="black", lw=1, label="Youden J = kept + rejected - 1")
    ac.plot(jo["t_opt"], jo["J_max"], "v", color="black", ms=6, label=f"max J = {jo['J_max']:.3f} at t = {jo['t_opt']:.3g}")
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
    for rule, col in (("flb_best", COLORS["sbw"]), ("d_dlr_min", COLORS["dlr"])):   # frozen DP0.2 J optima
        t = DP02_OPTIMA[rule]
        k = (~RULES[rule](df, t)).mean(); r = RULES[rule](rnd, t).mean()
        ad.plot(k, r, "D", color="none", ms=6, mec=col, mew=1.2,
                label=f"DP0.2 J optimum {'flb_best' if rule == 'flb_best' else '$d_{DLR}$'} = {t:g}: "
                      f"kept {100 * k:.1f}%, rejected {100 * r:.1f}%")
    ad.set_xlabel("fraction of real SNe kept (have a host)"); ad.set_ylabel("fraction of random positions rejected")
    kept_min = min((~RULES[r](df, t)).mean() for r, t in (("flb_best", cut), ("d_dlr_min", dmax)))
    ad.set_xlim(0.8 if kept_min > 0.85 else 0.5, 1.002); ad.set_ylim(0, 1.02)
    ad.legend(fontsize=7, loc="lower left"); ad.set_title("(d) hostless indicators compared", fontsize=9, loc="left")
    fig.suptitle(f"Hostless indicators on {LABEL['tier']}: {len(df):,} real SNe vs {len(rnd):,} random positions\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


# ---------------------------------------------------------------- credible sets
T_GRID = np.geomspace(1e-3, 0.5, 61)
T_MARK = (0.01, 0.05)
TOP_K = (1, 2, 3)


def credible_sets(df, grid=T_GRID):
    """Coverage of the true host by {surviving candidates with P_SBW > t} and the mean size of that set,
    over all SNe and over SNe whose true host survived the object-table cuts; plus the top-k sets."""
    tops = df["p_top"].to_numpy()
    L = max(1, max(len(x) for x in tops))
    P = np.zeros((len(tops), L), np.float32)
    for i, x in enumerate(tops):
        P[i, :len(x)] = x
    pt = p_true(df)
    kept = df["true_host_status"].to_numpy() == "kept"
    rows = []
    for t in grid:
        size = (P > t).sum(1)
        cov = pt > t
        rows.append({"t": float(t), "coverage": float(cov.mean()), "coverage_kept": float(cov[kept].mean()),
                     "mean_size": float(size.mean()), "frac_empty": float((size == 0).mean())})
    r = pd.DataFrame(rows)
    rank = df["rank_true"].to_numpy()
    nused = df["n_cand_used"].to_numpy()
    topk = [{"k": k, "coverage": float(((rank >= 1) & (rank <= k) & (pt > 0)).mean()),
             "coverage_kept": float(((rank >= 1) & (rank <= k) & (pt > 0))[kept].mean()),
             "mean_size": float(np.minimum(k, nused).mean())} for k in TOP_K]
    return r, topk


def knee(r):
    """Kneedle point of coverage versus mean set size: the grid point farthest above the chord joining the
    curve's ends, both axes scaled to [0, 1]."""
    x, y = r.mean_size.to_numpy(), r.coverage.to_numpy()
    xs, ys = (x - x.min()) / np.ptp(x), (y - y.min()) / np.ptp(y)
    d = ys - (ys[0] + (ys[-1] - ys[0]) * (xs - xs[0]) / (xs[-1] - xs[0]))
    return r.iloc[int(np.argmax(np.abs(d)))].to_dict()


def credible_set_metrics(df, out):
    r, topk = credible_sets(df)
    at = {f"t_{t:g}": credible_sets(df, [t])[0].iloc[0].to_dict() for t in T_MARK + (DP02_OPTIMA["credible_set_t"],)}
    kn = knee(r)
    res = {"n": len(df), "frac_true_host_kept": float((df.true_host_status == "kept").mean()),
           "at_marked_thresholds": at, "top_k": topk, "knee": kn,
           "dp02_knee_transfer": {"t": DP02_OPTIMA["credible_set_t"], **at[f"t_{DP02_OPTIMA['credible_set_t']:g}"],
                                  "own_knee_t": kn["t"],
                                  "coverage_minus_own_knee": at[f"t_{DP02_OPTIMA['credible_set_t']:g}"]["coverage"] - kn["coverage"]},
           "knee_definition": "grid point of coverage vs mean set size farthest from the chord between its ends (axes scaled to [0, 1])",
           "curve": r.round(6).to_dict(orient="list"), "settings": df.attrs}
    with open(out, "w") as f:
        json.dump(res, f, indent=2, default=float)


def credible_set(df, out):
    r, topk = credible_sets(df)
    kn = knee(r)
    ylo = min(0.8, np.floor(min(r.coverage.min(), min(k["coverage"] for k in topk)) * 20) / 20)
    fig, (aa, ab) = plt.subplots(1, 2, figsize=(8.6, 3.6))
    aa.plot(r.t, r.coverage, color=COLORS["sbw"], label="coverage, all SNe")
    aa.plot(r.t, r.coverage_kept, color=COLORS["sbw"], ls="--", lw=1, label="coverage, true host survived the cuts")
    aa.set_xscale("log"); aa.set_ylim(ylo, 1.002)
    aa.set_xlabel("threshold t (set = candidates with $P_{\\mathrm{SBW}}$ > t)"); aa.set_ylabel("fraction of SNe with the true host in the set")
    at = aa.twinx()
    at.plot(r.t, r.mean_size, color=COLORS["agree"], ls=":", label="mean set size (right axis)")
    at.set_ylabel("mean set size"); at.grid(False); at.spines["right"].set_visible(True)
    for t in T_MARK:
        aa.axvline(t, color=COLORS["ref"], ls=":", lw=1)
    h1, l1 = aa.get_legend_handles_labels(); h2, l2 = at.get_legend_handles_labels()
    aa.legend(h1 + h2, l1 + l2, fontsize=7, loc="lower left")
    aa.set_title("(a) coverage and set size versus threshold", fontsize=9, loc="left")
    ab.plot(r.mean_size, r.coverage, color=COLORS["sbw"], label="{P > t}, t = 0.001 ... 0.5")
    for t, mk in zip(T_MARK, ("o", "s")):
        q = credible_sets(df, [t])[0].iloc[0]
        ab.plot(q.mean_size, q.coverage, mk, color=COLORS["sbw"], mec="black", mew=0.6, ms=6,
                label=f"t = {t:g}: coverage {100 * q.coverage:.2f}%, size {q.mean_size:.2f}")
    ab.plot(kn["mean_size"], kn["coverage"], "v", color="black", ms=6, label=f"knee: t = {kn['t']:.3g}")
    ab.plot([k["mean_size"] for k in topk], [k["coverage"] for k in topk], "D-", color=COLORS["agree"], ms=4, lw=0.8,
            label="top-k, k = " + ", ".join(f"{k['k']}: {100 * k['coverage']:.2f}%" for k in topk))
    ab.set_xlabel("mean set size"); ab.set_ylabel("coverage"); ab.set_ylim(ylo, 1.002)
    ab.legend(fontsize=7, loc="lower right"); ab.set_title("(b) coverage versus set size", fontsize=9, loc="left")
    fig.suptitle(f"SBW credible sets on {LABEL['tier']} truth (true host not a surviving candidate = not covered)\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind", choices=["recovery", "calibration", "agreement_matrix", "agreement_vs_ddlr", "recovery_vs", "p_true_host",
                                    "hostless_indicator", "hostless_metrics", "credible_set", "credible_set_metrics"])
    ap.add_argument("--associations", required=True)
    ap.add_argument("--x", choices=sorted(XSPEC))
    ap.add_argument("--random", help="association table of random positions (hostless kinds)")
    ap.add_argument("--hostless-cut", choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--dlr-max", choices=sorted(DLR_MAX))
    ap.add_argument("--p-threshold", choices=sorted(P_TRUE))
    ap.add_argument("--tier", choices=sorted(TIER), default="dp02")
    ap.add_argument("--survivors-only", action="store_true",
                    help="agreement_matrix: real-SN panel on SNe whose true host passed the cat1 and cat2 cuts")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    LABEL["tier"] = TIER[a.tier]
    df = load(a.associations)
    if a.kind == "recovery":
        recovery(df, a.output, a.hostless_cut, load(a.random), a.p_threshold)
    elif a.kind == "calibration":
        calibration(df, a.output, a.hostless_cut)
    elif a.kind == "agreement_matrix":
        agreement_matrix(df, a.output, a.hostless_cut, load(a.random), a.p_threshold, survivors_only=a.survivors_only)
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
    elif a.kind == "credible_set":
        credible_set(df, a.output)
    elif a.kind == "credible_set_metrics":
        credible_set_metrics(df, a.output)


if __name__ == "__main__":
    main()
