"""Four-category sequential miss accounting for SBW and DLR: true host absent from the catalogue (cat1),
removed by object-table cuts (cat2), removed by the no-host cut (cat3), genuine association failure (cat4),
or success, each category counted only on the SNe that passed the previous one; for the headline
association and recomputed for every candidate_quality level. ``figure`` draws it from the json."""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

from dp02_fig import load
from dp02_lib import (COLORS, DLR_MAX, HOSTLESS_CUT, P_TRUE, PROFILES, QUALITY, RADIUS, associate_all,
                      load_candidates, miss_categories, wilson)

STAGES = ("absent", "cut", "hostless", "genuine")   # cat1..cat4, in the order they are applied
CATALOGUE = {"dp02": "DP0.2", "ou24": "OU2024"}


def frac(k, n):
    p, lo, hi = wilson(k, n)
    return {"n": int(k), "frac": float(p), "ci68": [float(lo), float(hi)]}


def account(df, cut, t, n_unmatched):
    """Counts and fractions (Wilson 68%) per method, of the full hosted sample and conditional on reaching
    the stage; the SNe the sample cut dropped because their true host has no catalogue object under
    truth_match count as cat1 (sub-category host_unmatched)."""
    mc = miss_categories(df, cut, t)
    n = len(df) + n_unmatched
    res = {"n": n, "n_associated": len(df), "n_host_unmatched": n_unmatched}
    for m in ("sbw", "dlr"):
        cat, sub = mc[m], mc[f"{m}_sub"]
        r, entered = {}, n
        for c in STAGES + ("success",):
            k = int((cat == c).sum()) + (n_unmatched if c == "absent" else 0)
            r[c] = {**frac(k, n), "n_entered": entered, "conditional": frac(k, entered)}
            entered -= k
        r["sub"] = {"absent": {"host_unmatched": n_unmatched, "not_primary_in_radius": int((sub == "absent").sum())},
                    "cut": {k: int((sub == k).sum()) for k in ("point_source", "unusable", "quality")},
                    "hostless": {"no_host_cut": int((sub == "hostless").sum())},
                    "genuine": {"p_below_threshold" if m == "sbw" else "other_closer":
                                int((sub == ("p_below_threshold" if m == "sbw" else "other_closer")).sum())}}
        if m == "dlr":
            r["sub"]["cut"]["unusable_note"] = "includes kept true hosts with no DLR ellipse"
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
    cat = CATALOGUE[a.tier]
    res = {"definition": {
        "order": "sequential: each category counts only SNe that passed every previous one; 'frac' is of all "
                 "hosted SNe, 'conditional' of the SNe that entered the stage ('n_entered')",
        "cat1_absent": f"true host has no {cat} object under truth_match (dropped by the sample cut, counted here) "
                       "or is not a primary object within the search radius",
        "cat2_cut": "true host removed before association: point source under extendedness_cut, no usable radius "
                    "or flux for the profile (always dropped), or failed candidate_quality; for DLR also a kept "
                    "true host with no DLR ellipse",
        "cat3_hostless_sbw": f"true host survived; SBW declares the SN hostless (no SBW host or flb_best < {cut:g})",
        "cat3_hostless_dlr": "true host survived; no candidate has d_DLR <= dlr_max",
        "cat4_genuine_sbw": f"a host is assigned but P_SBW(true) < {t:g}",
        "cat4_genuine_dlr": "DLR assigns another candidate with smaller d_DLR",
        "success_sbw": f"P_SBW(true) >= {t:g} and not SBW hostless", "success_dlr": "DLR host = true host"},
        "tier": a.tier,
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


STACK = [("success", "#4d9221"), ("genuine", "#d6604d"), ("hostless", "#92c5de"), ("cut", "#f4a582"),
         ("absent", "#bababa")]
SHORT = {"success": "success", "genuine": "cat4: genuine failure", "hostless": "cat3: no-host cut",
         "cut": "cat2: object-table cut", "absent": "cat1: not in catalogue"}


def figure(a):
    res = json.load(open(a.accounting))
    h = res["headline"]
    cat = CATALOGUE[res.get("tier", "dp02")]
    fig, (ax, ac, aq) = plt.subplots(1, 3, figsize=(13.2, 4.2), width_ratios=[1.25, 1, 1])
    lo = max(0.0, np.floor((min(h[m]["success"]["frac"] for m in ("sbw", "dlr")) - 0.03) * 20) / 20)
    for i, m in enumerate(("sbw", "dlr")):
        bottom = 0.0
        for c, color in STACK:
            f = h[m][c]["frac"]
            ax.bar(i, f, bottom=bottom, color=color, width=0.6, label=SHORT[c] if i == 0 else None)
            if bottom + f > lo and f > 0.002:
                ax.text(i + 0.33, max(bottom + f / 2, lo + 0.004), f"{100 * f:.2f}%", fontsize=7, va="center")
            bottom += f
    ax.set_xticks([0, 1], ["SBW\n($P_{\\mathrm{SBW}}$(true) $\\geq$ " + f"{res['p_true_threshold']:g})", "DLR"])
    ax.set_xlim(-0.5, 2.9); ax.set_ylim(lo, 1.0); ax.set_ylabel("fraction of hosted SNe")
    ax.legend(fontsize=7, loc="center right")
    ax.grid(axis="x", visible=False)
    ax.set_title(f"(a) outcome per method (N = {h['n']:,}, y axis from {lo:g})", fontsize=9, loc="left")
    # (b) conditional loss at each stage: of the SNe that reached it, the fraction lost there
    x = np.arange(len(STAGES))
    for m, color, dx in (("sbw", COLORS["sbw"], -0.18), ("dlr", COLORS["dlr"], 0.18)):
        f = np.array([h[m][c]["conditional"]["frac"] for c in STAGES])
        e = np.array([h[m][c]["conditional"]["ci68"] for c in STAGES])
        ac.bar(x + dx, f, width=0.36, color=color, label=m.upper())
        ac.errorbar(x + dx, f, yerr=[f - e[:, 0], e[:, 1] - f], fmt="none", ecolor="black", capsize=2, lw=0.8)
        for xi, fi in zip(x + dx, f):
            ac.text(xi, fi, f"{100 * fi:.1f}", ha="center", va="bottom", fontsize=6)
    ac.set_xticks(x, ["cat1\nnot in\ncatalogue", "cat2\nobject-table\ncut", "cat3\nno-host\ncut",
                      "cat4\ngenuine\nfailure"], fontsize=7)
    ac.set_ylabel("fraction lost, of the SNe reaching the stage"); ac.legend(fontsize=7)
    ac.grid(axis="x", visible=False)
    ac.set_title("(b) conditional loss per stage [%]", fontsize=9, loc="left")
    qs = list(res["quality_sweep"])
    xq = np.arange(len(qs))
    for m, color, dx in (("sbw", COLORS["sbw"], -0.06), ("dlr", COLORS["dlr"], 0.06)):
        for c, ls, mk in (("cut", "--", "s"), ("hostless", ":", "^"), ("genuine", "-", "o")):
            f = np.array([res["quality_sweep"][q][m][c]["frac"] for q in qs])
            e = np.array([res["quality_sweep"][q][m][c]["ci68"] for q in qs])
            aq.errorbar(xq + dx, f, yerr=[f - e[:, 0], e[:, 1] - f], color=color, ls=ls, marker=mk, ms=4, capsize=2,
                        lw=1.2, label=f"{m.upper()} {SHORT[c]}")
    aq.set_xticks(xq, qs); aq.set_xlabel("candidate_quality"); aq.set_ylabel("fraction of hosted SNe")
    aq.set_ylim(0, None); aq.legend(fontsize=6)
    aq.set_title("(c) cat2-cat4 as the object-table cuts tighten", fontsize=9, loc="left")
    s = res["settings"]
    fig.suptitle(f"Where the true host is lost on {cat}: catalogue, object-table cuts, no-host cut, or the association\n"
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
    c.add_argument("--tier", choices=sorted(CATALOGUE), default="dp02")
    c.add_argument("--workers", type=int, default=None)
    c.add_argument("--output", required=True)
    f = sub.add_parser("figure")
    f.add_argument("--accounting", required=True)
    f.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    compute(a) if a.mode == "compute" else figure(a)


if __name__ == "__main__":
    main()
