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

from dp02_lib import COLORS, calibration_curve, wilson

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5,
                     "legend.frameon": False, "figure.dpi": 120})


def load(path):
    t = pq.read_table(path)
    df = t.to_pandas()
    df.attrs = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if not k.startswith(b"pandas")}
    return df


def subtitle(df):
    a = df.attrs
    return f"{a.get('profile', '?')} | band {a.get('band', '?')} | {a.get('search_radius', '?')} | d_DLR<={a.get('dlr_max', '?')[1:]} ({a.get('dlr_radius', '?')}) | N={len(df)}"


def recovery(df, out):
    n = len(df)
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


def agreement_matrix(df, out):
    rows = [("SBW host = true", df.sbw_correct), ("SBW host wrong", ~df.sbw_correct)]
    cols = [("DLR host = true", df.dlr_correct), ("DLR host wrong", ~df.dlr_correct & ~df.dlr_hostless),
            ("DLR: no host", df.dlr_hostless)]
    n = len(df)
    M = np.array([[int((r & c).sum()) for _, c in cols] for _, r in rows])
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
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
    ax.set_xticks(range(3), [c for c, _ in cols]); ax.set_yticks(range(2), [r for r, _ in rows])
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind", choices=["recovery", "calibration", "agreement_matrix", "agreement_vs_ddlr", "recovery_vs", "p_true_host"])
    ap.add_argument("--associations", required=True)
    ap.add_argument("--x", choices=sorted(XSPEC))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = load(a.associations)
    if a.kind == "recovery":
        recovery(df, a.output)
    elif a.kind == "calibration":
        calibration(df, a.output)
    elif a.kind == "agreement_matrix":
        agreement_matrix(df, a.output)
    elif a.kind == "agreement_vs_ddlr":
        agreement_vs_ddlr(df, a.output)
    elif a.kind == "recovery_vs":
        recovery_vs(df, a.output, a.x)
    elif a.kind == "p_true_host":
        p_true_host(df, a.output)


if __name__ == "__main__":
    main()
