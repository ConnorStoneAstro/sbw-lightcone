"""Three-way comparison SBW / DLR / DELIGHT on the DP2 deep-drilling sample, drawn from
dp2_cross_method_metrics: for each pair of methods a 2x2 of hosted / hostless with the hosted-hosted cell
split into same host / different host, and the P_SBW bands of the DELIGHT host."""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dp02_fig import plt as _plt  # noqa: F401  (rcParams)
from dp02_lib import COLORS

PAIRS = ["SBW|DLR", "SBW|DELIGHT", "DLR|DELIGHT"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metrics", required=True, help="dp2_cross_method_metrics json")
    ap.add_argument("--subset", default="all", help="subset key (all, Legacy, PS1)")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    m = json.load(open(a.metrics))
    sub = m["subsets"][a.subset]
    n = sub["n"]
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.6), width_ratios=[1, 1, 1, 0.9])
    for ax, pair in zip(axes[:3], PAIRS):
        A, B = pair.split("|")
        c = sub["pairs"][pair]
        same, diff = c["same"]["k"], c["different"]["k"]
        a_only, b_only, both = c["a_only_hostless"]["k"], c["b_only_hostless"]["k"], c["both_hostless"]["k"]
        M = np.array([[same + diff, b_only], [a_only, both]], float)   # rows: A hosted / hostless; cols: B hosted / hostless
        ax.imshow(M / n, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        agree = "in SBW set" if A == "SBW" else "same host"
        txt = [[f"{agree} {same:,} ({100 * same / n:.1f}%)\n{'excluded' if A == 'SBW' else 'different'} {diff:,} ({100 * diff / n:.1f}%)", f"{b_only:,}\n{100 * b_only / n:.1f}%"],
               [f"{a_only:,}\n{100 * a_only / n:.1f}%", f"{both:,}\n{100 * both / n:.1f}%"]]
        for i in range(2):
            for j in range(2):
                ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=6.5, color="white" if M[i, j] / n > 0.5 else "black")
        ax.set_xticks([0, 1], [f"{B} hosted", f"{B} hostless"]); ax.set_yticks([0, 1], [f"{A} hosted", f"{A} hostless"])
        ax.grid(False)
        given = same / (same + diff) if same + diff else np.nan
        ax.set_title(f"{A} vs {B}: {'in SBW credible set' if A == 'SBW' else 'same host'} given both hosted {100 * given:.1f}%", fontsize=8)
    ax = axes[3]
    p = sub["p_sbw_of_delight_host"]
    keys = [k for k in ("gt_0p5", "t_to_0p5", "lt_t", "sbw_hostless", "not_a_candidate") if k in p]
    vals = [p[k]["frac"] for k in keys]
    ax.bar(range(len(keys)), vals, color=[COLORS["sbw"], COLORS["sbw"], COLORS["dlr"], COLORS["agree"], "#bababa"][:len(keys)])
    for i, v in enumerate(vals):
        ax.text(i, v, f"{100 * v:.1f}%", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(keys)), [k.replace("_", " ") for k in keys], fontsize=7, rotation=20)
    ax.set_ylabel("fraction of transients with a matched DELIGHT host"); ax.set_ylim(0, 1.05)
    ax.set_title(f"$P_{{\\mathrm{{SBW}}}}$ of the DELIGHT host (N = {p['n']:,})", fontsize=8)
    fig.suptitle(f"SBW, DLR and DELIGHT on the DP2 deep-drilling transients ({a.subset}, N = {n:,})", fontsize=9)
    fig.tight_layout(); fig.savefig(a.output); plt.close(fig)


if __name__ == "__main__":
    main()
