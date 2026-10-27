"""Three-category miss accounting for SBW and DLR: true host absent from the catalogue (cat1), removed by
object-table cuts (cat2), genuine association failure (cat3), or success; for the headline association and
recomputed for every candidate_quality level. ``figure`` draws it from the json."""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

from dp02_fig import load
from dp02_lib import (COLORS, DLR_MAX, HOSTLESS_CUT, MISS, P_TRUE, PROFILES, QUALITY, RADIUS,
                      associate_all, load_candidates, miss_categories, wilson)

SUBS = {"absent": ["host_unmatched", "not_primary_in_radius"],
        "cut": ["point_source", "unusable", "quality"],
        "genuine_sbw": ["p_below_threshold", "sbw_hostless"],
        "genuine_dlr": ["other_closer", "beyond_dlr_max"]}


def account(df, cut, t, n_unmatched):
    """Counts and fractions (Wilson 68%) per method; the SNe the sample cut dropped because their true host
    has no DP0.2 object under truth_match count as cat1 (sub-category host_unmatched)."""
    mc = miss_categories(df, cut, t)
    n = len(df) + n_unmatched
    res = {"n": n, "n_associated": len(df), "n_host_unmatched": n_unmatched}
    for m in ("sbw", "dlr"):
        cat, sub = mc[m], mc[f"{m}_sub"]
        r = {}
        for c in MISS:
            k = int((cat == c).sum()) + (n_unmatched if c == "absent" else 0)
            p, lo, hi = wilson(k, n)
            r[c] = {"n": k, "frac": float(p), "ci68": [float(lo), float(hi)]}
        subs = {"absent": {"host_unmatched": n_unmatched, "not_primary_in_radius": int((sub == "absent").sum())},
                "cut": {k: int((sub == k).sum()) for k in ("point_source", "unusable", "quality")}}
        if m == "dlr":
            subs["cut"]["unusable"] += int((sub == "unusable").sum())   # true host kept but no DLR ellipse
            subs["genuine"] = {k: int((sub == k).sum()) for k in SUBS["genuine_dlr"]}
        else:
            subs["genuine"] = {k: int((sub == k).sum()) for k in SUBS["genuine_sbw"]}
        r["sub"] = subs
        res[m] = r
    return res


def compute(a):
    cut, t = HOSTLESS_CUT[a.hostless_cut], P_TRUE[a.p_threshold]
    meta = pq.read_schema(a.sample).metadata
    precut = json.loads(meta[b"truth_match_precut"])
    sample = pq.read_table(a.sample).to_pandas(ignore_metadata=True)
    if "host_unmatched" in sample and sample.host_unmatched.any():   # cuts=all keeps them in the sample
        sample = sample[~sample.hostless & ~sample.host_missing]
        n_unmatched = 0
    else:
        n_unmatched = int(precut["n"] - len(sample))
    head = load(a.associations)
    res = {"definition": {
        "cat1_absent": "true host has no DP0.2 object under truth_match (dropped by the sample cut, counted here) "
                       "or is not a primary object within the search radius",
        "cat2_cut": "true host removed before association: point source under extendedness_cut, no usable radius "
                    "or flux for the profile (always dropped), or failed candidate_quality",
        "cat3_genuine_sbw": f"true host survived; P_SBW(true) < {t:g} or SBW hostless (flb_best < {cut:g})",
        "cat3_genuine_dlr": "true host survived; DLR assigns another candidate or none (d_DLR(true) > dlr_max)",
        "success_sbw": f"P_SBW(true) >= {t:g} and not SBW hostless", "success_dlr": "DLR host = true host"},
        "p_true_threshold": t, "sbw_hostless_cut": cut, "settings": head.attrs,
        "headline": account(head, cut, t, n_unmatched)}
    cand = load_candidates(a.candidates, bands=a.band)
    sweep = {}
    for q in QUALITY:
        df = associate_all(sample, cand, band=a.band, profile=a.profile, search_radius=a.search_radius,
                           dlr_max=a.dlr_max, dlr_radius=a.dlr_radius, extendedness=a.extendedness, quality=q,
                           workers=a.workers)
        r = account(df, cut, t, n_unmatched)
        r["frac_candidates_cut"] = {k: v for k, v in df.attrs.items() if k.startswith("frac_candidates_cut")}
        r["mean_p_sbw_true"] = float(np.nan_to_num(df.p_sbw_true.to_numpy(dtype=float)).mean())
        sweep[q] = r
    res["quality_sweep"] = sweep
    with open(a.output, "w") as f:
        json.dump(res, f, indent=2)


STACK = [("success", "#4d9221"), ("genuine", "#d6604d"), ("cut", "#f4a582"), ("absent", "#bababa")]
SHORT = {"success": "success", "genuine": "cat3: genuine failure", "cut": "cat2: object-table cut",
         "absent": "cat1: not in catalogue"}


def figure(a):
    res = json.load(open(a.accounting))
    h = res["headline"]
    fig, (ax, aq) = plt.subplots(1, 2, figsize=(9.6, 4.0), width_ratios=[1.3, 1])
    lo = max(0.0, np.floor((min(h[m]["success"]["frac"] for m in ("sbw", "dlr")) - 0.03) * 20) / 20)
    for i, m in enumerate(("sbw", "dlr")):
        bottom = 0.0
        for c, color in STACK:
            f = h[m][c]["frac"]
            ax.bar(i, f, bottom=bottom, color=color, width=0.6, label=SHORT[c] if i == 0 else None)
            mid = bottom + f / 2
            if bottom + f > lo and f > 0.002:
                ax.text(i + 0.33, max(mid, lo + 0.004), f"{100 * f:.2f}%", fontsize=7, va="center")
            bottom += f
    ax.set_xticks([0, 1], ["SBW\n($P_{\\mathrm{SBW}}$(true) $\\geq$ " + f"{res['p_true_threshold']:g})", "DLR"])
    ax.set_xlim(-0.5, 2.9); ax.set_ylim(lo, 1.0); ax.set_ylabel("fraction of hosted SNe")
    ax.legend(fontsize=7, loc="center right")
    ax.grid(axis="x", visible=False)
    ax.set_title(f"(a) outcome per method (N = {h['n']:,}, y axis from {lo:g})", fontsize=9, loc="left")
    qs = list(res["quality_sweep"])
    x = np.arange(len(qs))
    for m, color, dx in (("sbw", COLORS["sbw"], -0.06), ("dlr", COLORS["dlr"], 0.06)):
        for c, ls, mk in (("cut", "--", "s"), ("genuine", "-", "o")):
            f = np.array([res["quality_sweep"][q][m][c]["frac"] for q in qs])
            e = np.array([res["quality_sweep"][q][m][c]["ci68"] for q in qs])
            aq.errorbar(x + dx, f, yerr=[f - e[:, 0], e[:, 1] - f], color=color, ls=ls, marker=mk, ms=4, capsize=2,
                        lw=1.2, label=f"{m.upper()} {'cat2 (cut)' if c == 'cut' else 'cat3 (genuine)'}")
    aq.set_xticks(x, qs); aq.set_xlabel("candidate_quality"); aq.set_ylabel("fraction of hosted SNe")
    aq.set_ylim(0, None); aq.legend(fontsize=7)
    aq.set_title("(b) object-table cuts versus genuine failures", fontsize=9, loc="left")
    s = res["settings"]
    fig.suptitle("Where the true host is lost: catalogue, object-table cuts, or the association\n"
                 f"{s.get('profile')} | band {s.get('band')} | {s.get('search_radius')} | d_DLR<={s.get('dlr_max', 'd?')[1:]} "
                 f"({s.get('dlr_radius')}) | quality {s.get('candidate_quality')} | flb < {res['sbw_hostless_cut']:g}", fontsize=8)
    fig.tight_layout(); fig.savefig(a.output); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode", required=True)
    c = sub.add_parser("compute")
    c.add_argument("--sample", required=True)
    c.add_argument("--associations", required=True)
    c.add_argument("--candidates", required=True)
    c.add_argument("--profile", required=True, choices=PROFILES)
    c.add_argument("--band", required=True, choices=list("ugrizy"))
    c.add_argument("--search-radius", required=True, choices=sorted(RADIUS))
    c.add_argument("--dlr-max", required=True, choices=sorted(DLR_MAX))
    c.add_argument("--dlr-radius", required=True, choices=["half_light", "moments"])
    c.add_argument("--extendedness", required=True, choices=["ext05", "none"])
    c.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    c.add_argument("--p-threshold", required=True, choices=sorted(P_TRUE))
    c.add_argument("--workers", type=int, default=None)
    c.add_argument("--output", required=True)
    f = sub.add_parser("figure")
    f.add_argument("--accounting", required=True)
    f.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    compute(a) if a.mode == "compute" else figure(a)


if __name__ == "__main__":
    main()
