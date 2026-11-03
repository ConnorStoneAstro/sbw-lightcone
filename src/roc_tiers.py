"""Hostless ROC of every tier on one figure, from the *_hostless_metrics json files (their ``curves`` block):
random positions rejected versus real transients kept for SBW (flb_best) and DLR (d_DLR_min), with the fixed
operating points (sbw_hostless_cut, dlr_max) and the DP0.2 J optimum of the same universe marked."""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dp02_fig import plt as _plt  # noqa: F401  (rcParams)
from dp02_lib import COLORS, DLR_MAX, HOSTLESS_CUT

TIERS = [("dp02", "DC2 / DP0.2"), ("ou24", "OpenUniverse2024"), ("des", "DES-SN5YR"), ("dp2", "Rubin DP2")]


def point(m, rule, t):
    """Real kept and random rejected at threshold t, interpolated on the metric's own grid."""
    c = m["curves"]
    grid = np.array(c["flb_grid"] if rule != "d_dlr_min" else c["d_dlr_grid"])
    k, r = np.array(c[rule]["real_kept"]), np.array(c[rule]["random_rejected"])
    return float(np.interp(t, grid, k)), float(np.interp(t, grid, r))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    for key, _ in TIERS:
        ap.add_argument(f"--{key}", required=True, help=f"{key}_hostless_metrics json")
    ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--dlr-max", required=True, choices=sorted(DLR_MAX))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    cut, dmax = HOSTLESS_CUT[a.hostless_cut], DLR_MAX[a.dlr_max]
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.6), sharey=True)
    for ax, (key, name) in zip(axes, TIERS):
        m = json.load(open(getattr(a, key)))
        c = m["curves"]
        ax.plot(c["flb_best"]["real_kept"], c["flb_best"]["random_rejected"], color=COLORS["sbw"], label="SBW, flb_best < t")
        ax.plot(c["d_dlr_min"]["real_kept"], c["d_dlr_min"]["random_rejected"], color=COLORS["dlr"], label="DLR, $d_{DLR,\\min}$ > t")
        for rule, t, col, lab in (("flb_best", cut, COLORS["sbw"], f"flb_best = {cut:g}"), ("d_dlr_min", dmax, COLORS["dlr"], f"$d_{{DLR}}$ = {dmax:g}")):
            k, r = point(m, rule, t)
            ax.plot(k, r, "o", color=col, ms=6, mec="black", mew=0.6, label=f"{lab}: kept {100 * k:.1f}%, rej. {100 * r:.1f}%")
        opt = m["dp02_optima_transfer"]["dp02_optima"]
        for rule, col, lab in (("flb_best", COLORS["sbw"], "flb_best"), ("d_dlr_min", COLORS["dlr"], "$d_{DLR}$")):
            t = opt[rule]
            k, r = point(m, rule, t)
            ax.plot(k, r, "D", color="none", ms=6, mec=col, mew=1.2, label=f"DP0.2 J optimum {lab} = {t:.3g}")
        ax.set_xlim(0.5, 1.002); ax.set_ylim(0, 1.02)
        ax.set_title(f"{name}  (N = {m['n_real']:,} real, {m['n_random']:,} random)", fontsize=8)
        ax.set_xlabel("fraction of real transients kept")
        ax.legend(fontsize=6, loc="lower left")
    axes[0].set_ylabel("fraction of random positions rejected")
    fig.suptitle("Hostless indicators on every tier: random positions rejected versus real transients kept", fontsize=9)
    fig.tight_layout(); fig.savefig(a.output); plt.close(fig)


if __name__ == "__main__":
    main()
