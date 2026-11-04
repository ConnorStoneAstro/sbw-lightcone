"""Four-category sequential miss accounting for SBW and DLR: true host absent from the catalogue (cat1),
removed by object-table cuts (cat2), removed by the no-host cut (cat3), genuine association failure (cat4),
or success, each category counted only on the SNe that passed the previous one; for the headline
association and recomputed for every candidate_quality level. Where the true host is not available to the
method (cat1, cat2) the pipeline still returns something, so those SNe are split further by what it did:
called the SN hostless, or assigned a (necessarily wrong) host. The denominator is the sample (rung 2):
every detectable hosted SN with a host row within the offset cut, host-unmatched SNe included; the SNe the
peak-magnitude cut removed are reported from the sample's ladder metadata, not drawn. ``figure`` draws
the json."""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

from dp02_fig import load
from dp02_lib import (COLORS, DLR_MAX, HOSTLESS_CUT, P_TRUE, PLACEMENT, PLACEMENT_LABEL, PROFILES, QUALITY, RADIUS,
                      associate_all, load_candidates, miss_categories, placement, sample_mask, sbw_hostless, wilson)

STAGES = ("absent", "cut", "hostless", "genuine")   # cat1..cat4, in the order they are applied
CATALOGUE = {"dp02": "DP0.2", "ou24": "OU2024"}


def frac(k, n):
    p, lo, hi = wilson(k, n) if n else (np.nan, np.nan, np.nan)
    return {"n": int(k), "frac": float(p), "ci68": [float(lo), float(hi)]}


def account(df, cut, t):
    """Counts and fractions (Wilson 68%) per method, of the sample (rung 2) and conditional on reaching the
    stage. For cat1 and cat2 the pipeline outcome is split: hostless call or wrong host assigned."""
    mc = miss_categories(df, cut, t)
    n = len(df)
    tr = df["host_objectId"].to_numpy()
    res = {"n": n, "n_host_unmatched": int((tr < 0).sum())}
    for m in ("sbw", "dlr"):
        cat, sub = mc[m], mc[f"{m}_sub"]
        hostless_call = sbw_hostless(df, cut) if m == "sbw" else df["dlr_objectId"].to_numpy() < 0
        r, entered = {}, n
        for c in STAGES:
            sel = cat == c
            k = int(sel.sum())
            r[c] = {**frac(k, n), "n_entered": entered, "conditional": frac(k, entered)}
            if c in ("absent", "cut"):   # the true host is unavailable: what did the pipeline do instead?
                kh = int((sel & hostless_call).sum())
                r[c]["outcome"] = {"hostless_call": {**frac(kh, k), "of_sample": float(kh / n)},
                                   "wrong_host": {**frac(k - kh, k), "of_sample": float((k - kh) / n)}}
            entered -= k
        k = int((cat == "success").sum())   # conditional: of the SNe that entered cat4 (a host was assigned)
        r["success"] = {**frac(k, n), "n_entered": r["genuine"]["n_entered"],
                        "conditional": frac(k, r["genuine"]["n_entered"])}
        absent = cat == "absent"
        r["sub"] = {"absent": {"host_unmatched": int((absent & (tr < 0)).sum()),
                               "not_primary_in_radius": int((absent & (tr >= 0)).sum())},
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
    ladder = json.loads(meta[b"sample_ladder"]) if b"sample_ladder" in meta else None
    sample = pq.read_table(a.sample).to_pandas(ignore_metadata=True)
    sample = sample[sample_mask(sample)]           # a no-op unless cuts=all kept the flagged rows
    head = load(a.associations)
    attrs = head.attrs
    head = head[sample_mask(head)]
    head.attrs = attrs
    cat = CATALOGUE[a.tier]
    res = {"definition": {
        "denominator": "the sample (rung 2): hosted SNe with a host row, detectable under transient_peak_mag, within the "
                       "dp02_sample_cuts offset cut; SNe whose true host has no catalogue object are in it",
        "order": "sequential: each category counts only SNe that passed every previous one; 'frac' is of the "
                 "sample, 'conditional' of the SNe that entered the stage ('n_entered')",
        "cat1_absent": f"true host has no {cat} object under truth_match (host_unmatched) or is not a primary "
                       "object within the search radius",
        "cat2_cut": "true host removed before association: point source under extendedness_cut, no usable radius "
                    "or flux for the profile (always dropped), or failed candidate_quality; for DLR also a kept "
                    "true host with no DLR ellipse",
        "outcome": "for cat1 and cat2, what the pipeline did instead: hostless_call (SBW: no SBW host or flb_best "
                   f"< {cut:g}; DLR: no candidate with d_DLR <= dlr_max) or wrong_host (a host was assigned, which "
                   "cannot be the true one); 'frac' of the category, 'of_sample' of the sample",
        "cat3_hostless_sbw": f"true host survived; SBW declares the SN hostless (no SBW host, or flb_best < {cut:g}; "
                             "the same rule as the hostless indicator)",
        "cat3_hostless_dlr": "true host survived; no candidate has d_DLR <= dlr_max",
        "cat4_genuine_sbw": f"a host is assigned but P_SBW(true) < {t:g}",
        "cat4_genuine_dlr": "DLR assigns another candidate with smaller d_DLR",
        "success_sbw": f"P_SBW(true) >= {t:g} and not SBW hostless", "success_dlr": "DLR host = true host"},
        "tier": a.tier, "ladder": ladder,
        "p_true_threshold": t, "sbw_hostless_cut": cut, "settings": attrs,
        "placement_definition": {
            "detectable": "the sample: hosted transients with a host row, detectable under transient_peak_mag, within the offset cut",
            "hosted": "the method assigned a host (SBW: an SBW host with flb_best >= cut; DLR: a candidate with d_DLR <= dlr_max)",
            "rejected": "hostless call (the complement)",
            "true_host": "hosted and the true host is recovered (SBW: P_SBW(true) >= p_true_threshold, the true host is "
                         "in the credible set; DLR: the assigned host is the true host)",
            "wrong_host_miss": "hosted, the true host was available to the method (in the object table, kept by the cuts, "
                               "with a DLR ellipse for DLR) and not recovered (SBW: P_SBW(true) < p_true_threshold; DLR: "
                               "another object was assigned)",
            "wrong_host_not_in_table": "hosted but the true host was never a candidate (absent from the object table or "
                                       "removed by the object-table cuts), so the assigned host is wrong by construction",
            "rejected_true_host_present": "hostless call although the true host was available (the hostless cut fired)",
            "rejected_true_host_not_in_table": "hostless call and the true host was not available",
            "contamination": "wrong_host_miss + wrong_host_not_in_table over the hosted SNe",
            "fractions": "'frac' of the detectable sample, 'of_hosted' of the hosted SNe; Wilson 68% intervals"},
        "placement": placement(head, cut, t),
        "headline": account(head, cut, t)}
    cand = load_candidates(a.candidates, bands=a.band)
    sweep = {}
    for q in QUALITY:
        df = associate_all(sample, cand, band=a.band, profile=a.profile, search_radius=a.search_radius,
                           dlr_max=a.dlr_max, dlr_radius=a.dlr_radius, extendedness=a.extendedness, quality=q,
                           workers=a.workers)
        r = account(df, cut, t)
        r["frac_candidates_cut"] = {k: v for k, v in df.attrs.items() if k.startswith("frac_candidates_cut")}
        r["mean_p_sbw_true"] = float(np.nan_to_num(df.p_sbw_true.to_numpy(dtype=float)).mean())
        sweep[q] = r
    res["quality_sweep"] = sweep
    with open(a.output, "w") as f:
        json.dump(res, f, indent=2)


STACK = [("true_host", "#4d9221"), ("wrong_host_miss", "#d6604d"), ("wrong_host_not_in_table", "#67001f"),
         ("rejected_true_host_present", "#92c5de"), ("rejected_true_host_not_in_table", "#bababa")]
SHORT = {"true_host": "hosted, true host recovered",
         "wrong_host_miss": "hosted, true host available but not recovered",
         "wrong_host_not_in_table": "hosted, true host not in object table (wrong by construction)",
         "rejected_true_host_present": "rejected (true host was available)",
         "rejected_true_host_not_in_table": "rejected (true host not in object table)"}


def figure(a):
    """One panel: for SBW and DLR, how the detectable SNe are placed (hosted with the true host, hosted with
    a wrong host by either cause, or rejected), as fractions of the sample; the hosted-SN fractions of the
    two wrong-host outcomes are written beside the bar, since those are the contamination a user meets."""
    res = json.load(open(a.accounting))
    pl = res["placement"]
    n = pl["n_detectable"]
    t = res["p_true_threshold"]
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    lo = max(0.0, np.floor((min(pl[m]["true_host"]["frac"] for m in ("sbw", "dlr")) - 0.03) * 20) / 20)
    gap = 1.9 if a.legend == "outside" else 1.0   # room for the labels before the next bar
    xs = [0.0, gap]
    for i, m in enumerate(("sbw", "dlr")):
        bottom = 0.0
        for c, color in STACK:
            f = pl[m][c]["frac"]
            ax.bar(xs[i], f, bottom=bottom, color=color, width=0.6, label=SHORT[c] if i == 0 else None)
            label = f"{100 * f:.2f}%"
            if "of_hosted" in pl[m][c] and c != "true_host":
                label += f" ({100 * pl[m][c]['of_hosted']:.2f}% of hosted)"
            if bottom + f > lo and f > 0.0015:
                ax.text(xs[i] + 0.33, max(bottom + f / 2, lo + 0.004), label, fontsize=6.5, va="center")
            bottom += f
    ax.set_xticks(xs, [f"{m.upper()}\nhosted {100 * pl[m]['frac_hosted']:.1f}%\ncontam. {100 * pl[m]['contamination']['of_hosted']:.2f}%"
                           for m in ("sbw", "dlr")], fontsize=7.5)
    ax.set_ylim(lo, 1.0); ax.set_ylabel(f"fraction of detectable SNe (N = {n:,})")
    if a.legend == "outside":
        ax.set_xlim(-0.5, gap + 2.0)
        ax.legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
        ax.text(1.02, 0.66, f"recovered: SBW $P$(true host) $\\geq$ {t:g}\n(in the credible set); DLR host = true host\n\n"
                            "contam. = not recovered, either cause,\nover the hosted SNe", transform=ax.transAxes, fontsize=6.5, va="top")
    else:
        ax.set_xlim(-0.5, 3.1)
        ax.legend(fontsize=6.5, loc="center right")
    ax.grid(axis="x", visible=False)
    ax.set_title("How SNe are placed with a host", fontsize=10, loc="left")
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
    f.add_argument("--legend", choices=["inside", "outside"], default="inside",
                   help="outside: legend and note to the right of the axes (paper layout)")
    f.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    compute(a) if a.mode == "compute" else figure(a)


if __name__ == "__main__":
    main()
