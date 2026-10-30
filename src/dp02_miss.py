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
from dp02_lib import (COLORS, DLR_MAX, HOSTLESS_CUT, P_TRUE, PROFILES, QUALITY, RADIUS, associate_all,
                      load_candidates, miss_categories, sample_mask, sbw_hostless, wilson)

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
        "cat3_hostless_sbw": f"true host survived; SBW declares the SN hostless (no SBW host or flb_best < {cut:g})",
        "cat3_hostless_dlr": "true host survived; no candidate has d_DLR <= dlr_max",
        "cat4_genuine_sbw": f"a host is assigned but P_SBW(true) < {t:g}",
        "cat4_genuine_dlr": "DLR assigns another candidate with smaller d_DLR",
        "success_sbw": f"P_SBW(true) >= {t:g} and not SBW hostless", "success_dlr": "DLR host = true host"},
        "tier": a.tier, "ladder": ladder,
        "p_true_threshold": t, "sbw_hostless_cut": cut, "settings": attrs,
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


STACK = [("success", "#4d9221"), ("genuine", "#d6604d"), ("hostless", "#92c5de"), ("cut", "#f4a582"),
         ("absent", "#bababa")]
SHORT = {"success": "success", "genuine": "cat4: genuine failure", "hostless": "cat3: no-host cut",
         "cut": "cat2: object-table cut", "absent": "cat1: not in catalogue"}
OUTCOME = [("wrong_host", "#67001f", "pipeline assigned a wrong host"), ("hostless_call", "#2166ac", "pipeline called it hostless")]


def figure(a):
    res = json.load(open(a.accounting))
    h = res["headline"]
    cat = CATALOGUE[res.get("tier", "dp02")]
    fig, (ax, ac) = plt.subplots(1, 2, figsize=(9.6, 4.4), width_ratios=[1.35, 1])
    lo = max(0.0, np.floor((min(h[m]["success"]["frac"] for m in ("sbw", "dlr")) - 0.03) * 20) / 20)
    W, WI = 0.6, 0.22
    for i, m in enumerate(("sbw", "dlr")):
        bottom = 0.0
        for c, color in STACK:
            f = h[m][c]["frac"]
            ax.bar(i, f, bottom=bottom, color=color, width=W, label=SHORT[c] if i == 0 else None)
            label = f"{100 * f:.2f}%"
            if "outcome" in h[m][c] and f > 0:
                # bar in bar: the pipeline's outcome for SNe whose true host it could not have found
                b = bottom
                for key, col, lab in OUTCOME:
                    fo = h[m][c]["outcome"][key]["of_sample"]
                    ax.bar(i, fo, bottom=b, color=col, width=WI, label=lab if (i == 0 and c == "absent") else None)
                    b += fo
                label += f"\nwrong host {100 * h[m][c]['outcome']['wrong_host']['of_sample']:.2f}%"
            if bottom + f > lo and f > 0.002:
                ax.text(i + 0.33, max(bottom + f / 2, lo + 0.004), label, fontsize=6.5, va="center")
            bottom += f
    ax.set_xticks([0, 1], ["SBW\n($P_{\\mathrm{SBW}}$(true) $\\geq$ " + f"{res['p_true_threshold']:g})", "DLR"])
    ax.set_xlim(-0.5, 2.9); ax.set_ylim(lo, 1.0); ax.set_ylabel("fraction of the sample (rung 2)")
    ax.legend(fontsize=6.5, loc="center right")
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
    ac.set_title("(b) conditional loss per stage (labels in %)", fontsize=9, loc="left")
    s = res["settings"]
    lad = res.get("ladder") or {}
    peak = f" | peak {lad.get('peak_mag_cut', '?')}" if lad else ""
    fig.suptitle(f"Where the true host is lost on {cat}: catalogue, object-table cuts, no-host cut, or the association\n"
                 f"{s.get('profile')} | band {s.get('band')} | {s.get('search_radius')} | d_DLR<={s.get('dlr_max', 'd?')[1:]} "
                 f"({s.get('dlr_radius')}) | quality {s.get('candidate_quality')} | flb < {res['sbw_hostless_cut']:g}{peak}", fontsize=8)
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
