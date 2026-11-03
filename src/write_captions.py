"""Write a long-form caption (markdown) beside every paper-level figure and table, from the result files:
what the figure shows, which outputs and universe it comes from (paths), the sample it stands on with its N,
the decisions behind it (from the manifests and the spec), and the numbers the text should quote. Run from
the project root after a materialization: ``python src/write_captions.py`` writes ``captions/*.md``. It is a
project tool, not an output: nothing here is a recipe, and it reads results it does not make."""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "captions"
POINTER = 300   # bytes: a git-annex pointer, the content is elsewhere

FIXED = ("Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); "
         "the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with "
         "d_DLR <= 4 on the second-moment ellipse, else DLR hostless.")
LADDER = ("Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, "
          "with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within "
          "3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); "
          "the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted "
          "transient holds the *true host*, a *wrong host because the true host was not in the object table* "
          "(absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the "
          "true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host "
          "survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.")


def readable(path):
    p = ROOT / path
    return p.exists() and p.stat().st_size > POINTER


def load_json(path):
    return json.load(open(ROOT / path)) if readable(path) else None


def manifest(universe, oid):
    p = ROOT / "results" / universe / f".{oid}.manifest.json"
    return json.load(open(p)) if p.exists() else {}


def pct(x, d=2):
    return "n/a" if x is None or x != x else f"{100 * x:.{d}f} %"


def ci(v, d=2):
    return "n/a" if not v else f"[{100 * v[0]:.{d}f}, {100 * v[1]:.{d}f}]"


def head(title, outputs, sample, decisions_note):
    lines = [f"# {title}", ""]
    lines.append("## Outputs")
    for universe, oid, fmt in outputs:
        m = manifest(universe, oid)
        path = f"results/{universe}/{oid}.{fmt}"
        where = "content on this machine" if readable(path) else "annex pointer here; content on the machine that made it"
        lines.append(f"- `{path}` ({where}); made at commit `{m.get('git_sha', '?')[:8]}` with `{m.get('recipe', '?')}`")
        if m.get("decisions"):
            lines.append("  - decisions in the recipe: " + ", ".join(f"{k} = {v}" for k, v in m["decisions"].items()))
    lines += ["", "## Sample", "", sample, "", LADDER, "", "## Decisions behind it", "", decisions_note, "", FIXED, ""]
    return lines


def settings_line(j):
    s = (j or {}).get("settings") or {}
    keys = ("profile", "band", "search_radius", "dlr_max", "dlr_radius", "extendedness", "candidate_quality",
            "truth_match", "dp02_sample_cuts", "transient_peak_mag", "dp2_sample", "dp2_candidate_catalogue",
            "reference_match", "exclude_suspect")
    return ", ".join(f"{k} = {s[k]}" for k in keys if k in s) or "(settings block missing)"


# ---------------------------------------------------------------- figures
def fig1():
    lines = head("Fig 1. How SNe are placed with a host (DP0.2 and OU2024)",
                 [("baseline", "dp02_miss_accounting_fig", "pdf"), ("baseline", "ou24_miss_accounting_fig", "pdf"),
                  ("baseline", "dp02_miss_accounting", "json"), ("baseline", "ou24_miss_accounting", "json")],
                 "Every detectable SN of the simulation (the sample); this is the only figure that shows the whole sample, "
                 "so the reader sees the catalogue's losses and the method's losses side by side.",
                 "One stacked bar per method. From the bottom: hosted with the true host (green); hosted with a wrong host "
                 "although the true host was available (red: the method's real misses; for SBW the most probable host is "
                 "wrong, and the json counts how many of those still hold the true host at P >= 0.05); hosted with a wrong "
                 "host because the true host was not in the object table (dark red: a catalogue effect, wrong by "
                 "construction); rejected although the true host was available (light blue: the hostless cut fired); "
                 "rejected with the true host not in the object table (grey). Beside the bar each fraction of the sample "
                 "and, for the two wrong-host outcomes, the fraction of the *hosted* SNe, which is the contamination a "
                 "user of the hosted sample meets. Under the bar: hosted fraction and contamination. The transients "
                 "dropped as undetectable (peak r >= 24 or outside the imaged years) are not in the figure: their count "
                 "is the survey selection function and is quoted from the ladder below. The choice to keep host-unmatched "
                 "SNe in the sample (31 Oct 2026) is what makes the dark-red and grey segments possible: before, those SNe "
                 "were dropped before association and nobody knew what the pipeline did with them.")
    lines += ["## Numbers", ""]
    for tier, name in (("dp02", "DC2 / DP0.2"), ("ou24", "OpenUniverse2024")):
        j = load_json(f"results/baseline/{tier}_miss_accounting.json")
        if not j or "placement" not in j:
            lines.append(f"- {name}: result file not readable here"); continue
        lad, pl = j.get("ladder") or {}, j["placement"]
        lines.append(f"### {name} (settings: {settings_line(j)})")
        lines.append(f"- ladder: transients {lad.get('transients'):,}, hosted {lad.get('hosted'):,}, with host row {lad.get('host_row'):,}, "
                     f"detectable {lad.get('detectable'):,} (faint dropped {lad.get('faint_dropped'):,}), beyond the offset cut "
                     f"{lad.get('far_dropped'):,}; sample {lad.get('sample'):,} of which host-unmatched {lad.get('host_unmatched_in_sample'):,}")
        lines.append(f"- true host not in the object table: {pl['n_true_host_not_in_table']:,} ({pct(pl['n_true_host_not_in_table'] / pl['n_detectable'])}): "
                     f"absent {pl['n_true_host_absent']:,}, removed by the cuts {pl['n_true_host_cut']:,}")
        for m in ("sbw", "dlr"):
            r = pl[m]
            lines.append(f"- {m.upper()}: hosted {r['n_hosted']:,} ({pct(r['frac_hosted'])}), rejected {r['n_rejected']:,}; "
                         f"true host {pct(r['true_host']['frac'])} of sample; wrong host, true host available {r['wrong_host_miss']['n']:,} "
                         f"({pct(r['wrong_host_miss']['frac'])} of sample, {pct(r['wrong_host_miss']['of_hosted'])} of hosted"
                         + (f", of which {r['wrong_host_miss'].get('true_host_p_ge_t')} still hold the true host at P >= 0.05" if m == "sbw" else "")
                         + f"); wrong host, true host not in table {r['wrong_host_not_in_table']['n']:,} ({pct(r['wrong_host_not_in_table']['of_hosted'])} of hosted; "
                         f"absent {r['wrong_host_not_in_table']['absent']}, cut {r['wrong_host_not_in_table']['cut']}); rejected with true host available "
                         f"{r['rejected_true_host_present']['n']:,} ({pct(r['rejected_true_host_present']['frac'])}); rejected, true host not in table "
                         f"{r['rejected_true_host_not_in_table']['n']:,}; contamination {pct(r['contamination']['of_hosted'])} {ci(r['contamination']['ci68'])} of hosted")
        lines.append("")
    return lines


def fig2():
    lines = head("Fig 2. Candidate-level reliability of P_SBW (DP0.2 and OU2024; DES in text)",
                 [("baseline", "dp02_calibration", "pdf"), ("baseline", "ou24_calibration", "pdf"), ("baseline", "des_calibration", "pdf")],
                 "Hosted SNe (SBW assigned a host at flb_best >= 0.01) whose true host is usable (rung 4), so some candidate "
                 "is the right one; every candidate with P_SBW > 0.001 of those SNe enters, binned by its P_SBW. SNe that SBW "
                 "rejects and SNe whose true host was not in the object table are left out: SBW's probabilities are not used "
                 "on the first and cannot be right on the second.",
                 "Per bin, the fraction of candidates that are the true host against the bin's mean P_SBW; a calibrated "
                 "probability sits on the 1:1 line in every bin, not only for the most probable host. Log-log over the full "
                 "range (left), linear above 0.5 (right), candidate and true-host counts per bin below. The calibration error "
                 "is the probability mass per SN that lands in the wrong bin: sum over bins of |true hosts - summed P| / N. "
                 "The earlier top-host-only curve conditioned on the ranking and hid the low-P bins where SBW decides how "
                 "much mass faint neighbours get (30 Oct 2026). The true-host probability is rounded to float32 like the "
                 "candidate list so both histograms share the same bin membership.")
    lines += ["## Numbers", ""]
    for tier, name in (("dp02", "DC2 / DP0.2"), ("ou24", "OpenUniverse2024"), ("des", "DES-SN5YR (reference = published host)")):
        j = load_json(f"results/baseline/{tier}_recovery.json")
        if not j:
            lines.append(f"- {name}: result file not readable here"); continue
        b = j.get("p_true_bands", {})
        lines.append(f"- {name}: calibration error {j['expected_calibration_error']:.4f} per SN on {j['calibration_n_sne']:,} SNe; rung 4 N = {j['n']:,}; "
                     f"mean P_SBW(reference) {pct(j['mean_p_sbw_true'])} +- {100 * j.get('mean_p_sbw_true_se', 0):.2f}; "
                     f"P bands: > 0.5 {pct(b.get('p_gt_0p5'))}, 0.05..0.5 {pct(b.get('p_t_to_0p5'))}, < 0.05 {pct(b.get('p_lt_t'))}, SBW hostless {pct(b.get('sbw_hostless'))}")
    return lines


def fig3():
    lines = head("Fig 3. SBW probability of the true host against DLR, on SNe whose true host is usable (DP0.2 and OU2024)",
                 [("baseline", "dp02_agreement_matrix", "pdf"), ("baseline", "ou24_agreement_matrix", "pdf")],
                 "Rung 4: SNe whose true host passed the object-table cuts for both methods, so every cell is an association "
                 "outcome. The random-position panel beside it uses the ten-per-SN annulus draws under the same cuts.",
                 "Rows: P_SBW(true host) > 0.5, in [0.05, 0.5], < 0.05 (excluded), or SBW hostless at flb_best < 0.01; columns: "
                 "DLR host = true, DLR another object, DLR hostless. Counts and fractions of rung 4 in every cell. The "
                 "random 2x2 shows how often each method correctly rejects a position with no host. The 'both wrong, same "
                 "object' tail of DC2 comes from its untruncated Sersic placement; the physical offset cut (3 R_e) removes most "
                 "of it (dp02_offset_diagnostic).")
    lines += ["## Numbers", ""]
    for tier, name in (("dp02", "DC2 / DP0.2"), ("ou24", "OpenUniverse2024")):
        j = load_json(f"results/baseline/{tier}_recovery.json")
        if not j:
            lines.append(f"- {name}: result file not readable here"); continue
        m, r = j["agreement_matrix"], j["random_matrix"]
        n = m["n"]
        cells = ", ".join(f"{k} {v:,} ({pct(v / n)})" for k, v in m.items() if "/" in k)
        lines.append(f"- {name}: N = {n:,}; {cells}; SBW argmax = true {pct(j['sbw_argmax_recovery'])} {ci(j['sbw_argmax_recovery_ci68'])}, "
                     f"DLR = true {pct(j['dlr_recovery'])} {ci(j['dlr_recovery_ci68'])}, SBW most probable = DLR {pct(j['agreement'])}; "
                     f"random positions N = {r['n']:,}: SBW rejects {pct(r['sbw_rejected'])}, DLR rejects {pct(r['dlr_rejected'])}")
    return lines


def fig4():
    lines = head("Fig 4. Recovery at fixed host magnitude, DP0.2 against OU2024",
                 [("baseline", "ou24_dp02_comparison", "pdf"), ("baseline", "ou24_recovery_by_host_mag", "json")],
                 "Rung 4 of both simulations (true host usable), binned by the true host's r magnitude; a grey line per tier "
                 "gives the fraction of sample hosts that reach rung 4 in each bin (the catalogue's selection function).",
                 "OU2024 hosts are drawn about 1 dex too light (OU2024 paper section 3.2.4), so the population-averaged "
                 "numbers of the two simulations are not comparable; at fixed host magnitude the association is the same "
                 "geometry and photometry problem in both. Curves: SBW mean P_SBW(true) and DLR fraction correct, DP0.2 solid, "
                 "OU2024 dashed (all z, and z < 1 to match DC2's range); the third panel is the host-r distributions. The "
                 "text box gives OU2024 reweighted to DP0.2's host-r distribution and both tiers at host r < 24.")
    lines += ["## Numbers", ""]
    j = load_json("results/baseline/ou24_recovery_by_host_mag.json")
    if j and "rung4" in j:
        for r in ("rung4", "rung2"):
            for z in ("all_z", "z_lt_1"):
                a, b = j[r][z]["all"], j[r][z]["host_r_lt_24"]
                lines.append(f"- OU2024 {r} ({j[r]['label']}), {z}: N {a['n']:,}, mean P {pct(a['mean_p_sbw_true'])}, DLR {pct(a['dlr_recovery'])}; "
                             f"host r < 24: N {b['n']:,}, mean P {pct(b['mean_p_sbw_true'])}, DLR {pct(b['dlr_recovery'])}")
    else:
        lines.append("- ou24_recovery_by_host_mag not readable here")
    d = load_json("results/baseline/dp02_recovery.json")
    if d:
        lines.append(f"- DP0.2 rung 4: N {d['n']:,}, mean P {pct(d['mean_p_sbw_true'])}, DLR {pct(d['dlr_recovery'])}")
    return lines


def fig5():
    lines = head("Fig 5. Hostless indicators on every tier (ROC)",
                 [("baseline", "hostless_roc_tiers", "pdf")] + [("baseline", f"{t}_hostless_metrics", "json") for t in ("dp02", "ou24", "des", "dp2")],
                 "Real transients on rung 4 where a reference exists (DP0.2, OU2024, DES), all 690 rc_salt transients on DP2; "
                 "random positions are ten per transient in the 45-90 arcsec annulus around it on every tier (2 Nov 2026; before, "
                 "the designs differed by tier and so did the random baselines).",
                 "For SBW the hostless statistic is frac_light_beyond of the most probable host (flb_best): the fraction of "
                 "that host's light that lies farther from its centre than the transient; for DLR it is the smallest d_DLR. "
                 "Each curve sweeps the threshold; filled circles are the fixed operating points (flb 0.01, d_DLR 4) and open "
                 "diamonds the Youden-J optima of the DP0.2 tier of the same universe (read from dp02_hostless_metrics, not "
                 "frozen constants), so the other three tiers are held-out tests of whether a threshold tuned on DC2 transfers. "
                 "Random positions get a P > 0.5 most probable host on every tier almost always, so hostlessness must come from "
                 "flb, never from P.")
    lines += ["## Numbers", ""]
    for t, name in (("dp02", "DC2 / DP0.2"), ("ou24", "OpenUniverse2024"), ("des", "DES-SN5YR"), ("dp2", "Rubin DP2 rc_salt")):
        j = load_json(f"results/baseline/{t}_hostless_metrics.json")
        if not j:
            lines.append(f"- {name}: result file not readable here"); continue
        y = j["youden_J"]
        tr = j["dp02_optima_transfer"]
        pts = {p["threshold"]: p for p in tr["flb_best"]["at"]}
        dpt = {p["threshold"]: p for p in tr["d_dlr_min"]["at"]}
        p01 = pts.get(0.01) or next(iter(pts.values()))
        d4 = dpt.get(4.0) or next(iter(dpt.values()))
        lines.append(f"- {name}: N real {j['n_real']:,} (sample {j.get('n_sample', j['n_real']):,}), random {j['n_random']:,}; "
                     f"at flb 0.01 kept {pct(p01['real_kept'], 1)}, random rejected {pct(p01['random_rejected'], 1)}, J {p01['J']:.3f}; "
                     f"at d_DLR 4 kept {pct(d4['real_kept'], 1)}, rejected {pct(d4['random_rejected'], 1)}, J {d4['J']:.3f}; "
                     f"own J optimum flb {y['sbw_flb_best']['t_opt']:.3g} (J {y['sbw_flb_best']['J_max']:.3f}), d_DLR {y['dlr_d_dlr_min']['t_opt']:.3g} (J {y['dlr_d_dlr_min']['J_max']:.3f}); "
                     f"DP0.2 optima used: flb {tr['dp02_optima']['flb_best']:.3g}, d_DLR {tr['dp02_optima']['d_dlr_min']:.3g} (source {tr['dp02_optima'].get('source', '?')}); "
                     f"random positions landing on a P > 0.5 host {pct(j['random_p_sbw_gt_0p5'], 1)}")
    return lines


def fig6():
    lines = head("Fig 6. DES-SN5YR: SBW and DLR against the published host",
                 [("baseline", "des_agreement_matrix", "pdf"), ("baseline", "des_recovery", "json")],
                 "The 1,623 DES-SN5YR SNe Ia with a published host; the matrix is on rung 4 (published host matched to a DES "
                 "DR2 object within 1 arcsec, kept by the object-table cuts, with a DLR ellipse). The reference is itself a "
                 "DLR assignment on deeper SN-field coadds, so this is agreement, not recovery.",
                 "Same layout as Fig 3 with the published host as reference; candidates are DES DR2 objects within 30 arcsec "
                 "with SExtractor light models; random positions ten per SN in the annulus. HOSTGAL2 (the published second "
                 "host) explains most of the remaining 'other object' picks; the SNe whose published host is absent from DR2 "
                 "or classed as a point source are the catalogue part of the ladder and appear only in the numbers.")
    lines += ["## Numbers", ""]
    j = load_json("results/baseline/des_recovery.json")
    if not j:
        lines.append("- des_recovery not readable here"); return lines
    a, h2, st = j["agreement_with_published"], j["hostgal2"], j["published_host_status"]
    lines.append(f"- settings: {settings_line(j)}")
    lines.append(f"- ladder: {j['ladder']}; rung 4 N = {j['n']:,}")
    lines.append(f"- SBW most probable = published {pct(a['sbw_argmax']['frac'])} {ci(a['sbw_argmax']['ci68'])}; hosted and = published {pct(a['sbw_argmax_hosted']['frac'])}; "
                 f"DLR = published {pct(a['dlr']['frac'])} {ci(a['dlr']['ci68'])}; mean P_SBW(published) {pct(a['mean_p_sbw_published'])}")
    lines.append(f"- published host status over the sample: {st}; matched in DR2 {pct(a['published_matched_in_dr2']['frac'])}")
    lines.append(f"- HOSTGAL2: present {h2['n_present']}, matched {h2['n_matched']}; of SBW picks that differ from the published host, HOSTGAL2 "
                 f"{pct(h2['of_sbw_ne_published_hosted_is_host2']['frac'])}; DLR {pct(h2['of_dlr_ne_published_hosted_is_host2']['frac'])}")
    m, r = j["agreement_matrix"], j["random_matrix"]
    lines.append("- matrix cells: " + ", ".join(f"{k} {v:,}" for k, v in m.items() if "/" in k) + f"; random N {r['n']:,}, SBW rejects {pct(r['sbw_rejected'])}, DLR {pct(r['dlr_rejected'])}")
    lines.append(f"- calibration error {j['expected_calibration_error']:.4f} on {j['calibration_n_sne']:,} SNe")
    pb = j.get("per_band_profile", {})
    if "r/exponential" in pb:
        lines.append(f"- per band/profile (same denominator as the headline): r/exponential SBW argmax {pct(pb['r/exponential']['sbw_argmax_agrees'])}, DLR {pct(pb['r/exponential']['dlr_agrees'])}")
    return lines


def fig7():
    lines = head("Fig 7. Rubin DP2 without a reference: SBW against DLR",
                 [("baseline", "dp2_agreement_matrix", "pdf"), ("baseline", "dp2_disagreement_rate", "json")],
                 "The 690 rc_salt transients (SALT3-fit SNe Ia from the DP2 release candidate) with DP2 Object candidates "
                 "within 30 arcsec and MultiProFit light models; random positions ten per transient in the annulus.",
                 "Rows: P_SBW of the DLR host > 0.5 / in [0.05, 0.5] / < 0.05, or SBW hostless; columns: the DLR host is SBW's "
                 "most probable host / another object / DLR hostless. Without truth the question is whether the two methods "
                 "name the same object when both assign one, and where their hostless calls differ (compact faint hosts "
                 "whose moment ellipse is tiny).")
    lines += ["## Numbers", ""]
    j = load_json("results/baseline/dp2_disagreement_rate.json")
    if not j:
        lines.append("- dp2_disagreement_rate not readable here"); return lines
    h, r = j["headline"], j["random_positions"]
    lines.append(f"- settings: {settings_line(j)}")
    lines.append(f"- N {h['n']}: same host {pct(h['same_host'])}, different host {pct(h['different_host'])}, SBW hostless only {pct(h['sbw_hostless_only'])}, "
                 f"DLR hostless only {pct(h['dlr_hostless_only'])}, both hostless {pct(h['both_hostless'])}; SBW hostless {pct(h['sbw_hostless'])}, DLR hostless {pct(h['dlr_hostless'])}; "
                 f"same host given DLR hosted {pct(h['same_host'] / (1 - h['dlr_hostless']))}; mean P_SBW(DLR host) {pct(h['mean_p_sbw_of_dlr_host'])}; "
                 f"P(DLR host) < 0.05 {pct(h['frac_p_sbw_of_dlr_host_lt_t'])}"
                 + (f"; P of the most probable host: median {h['p_sbw_max_median_hosted']:.3f}, > 0.9 for {pct(h['p_sbw_max_gt_0p9_hosted'])}" if "p_sbw_max_median_hosted" in h else ""))
    lines.append(f"- random positions N {r['n']:,}: SBW hostless {pct(r['sbw_hostless'])}, DLR hostless {pct(r['dlr_hostless'])}")
    return lines


def fig8():
    lines = head("Fig 8. Band and light-model consistency of the SBW host on real data (DP2)",
                 [("baseline", "dp2_band_consistency", "pdf"), ("baseline", "dp2_disagreement_rate", "json")],
                 "The 690 rc_salt DP2 transients, the association recomputed in every band u-y and for every light model "
                 "(exponential fit, free-n Sersic fit, moment Gaussian) from dp2_association_sweep.",
                 "(a) per band, how often the SBW most probable host and the DLR host differ from the r-band ones when both "
                 "name one; (b) the SBW-DLR categories per band; (c) hostless fractions per band and light model; (d) SBW "
                 "host agreement between light models per band. Without truth this is the stability of the answer, not its "
                 "correctness; u is the exception for lack of flux.")
    lines += ["## Numbers", ""]
    j = load_json("results/baseline/dp2_disagreement_rate.json")
    if not j:
        lines.append("- dp2_disagreement_rate not readable here"); return lines
    prof = (j.get("settings") or {}).get("profile", "exponential")
    bf = j.get("band_flip", {}).get(prof, {})
    lines.append(f"- SBW host changes vs r band ({prof}), given hosted in both: " + ", ".join(f"{b} {pct(v['sbw_host_changes'])}" for b, v in bf.items()))
    lines.append("- hostless status changes vs r band: " + ", ".join(f"{b} {pct(v['sbw_hostless_status_changes'])}" for b, v in bf.items()))
    pa = j.get("profile_agreement", {})
    lines.append("- same SBW host between light models (r band, both hosted): " + ", ".join(f"{k} {pct(v['same_sbw_host'])} (N {v['n_both_hosted']})" for k, v in pa.items()))
    return lines


def fig9():
    lines = head("Fig 9. SBW, DLR and DELIGHT on the DP2 deep-drilling transients",
                 [("dp2_ddf_clean", "dp2_cross_method_matrix", "pdf"), ("dp2_ddf_clean", "dp2_cross_method_metrics", "json")],
                 "The cleaned DDF ParSNIP sample (dp2_sample = parsnip_ddf_clean: the DDF transients with a DELIGHT host "
                 "from Y.-J. Kang, minus those on a pre-existing Legacy point source that is a Gaia star or a WISE AGN); "
                 "DELIGHT's host is matched to a DP2 Object within 1 arcsec (reference_match = pos_1arcsec).",
                 "For each pair of methods a 2x2 of hosted / hostless with the hosted-hosted cell split into same host / "
                 "different host; the title gives same host given both hosted. The last panel is the P_SBW band of the "
                 "DELIGHT host over the transients whose DELIGHT position matched a candidate. DELIGHT has no hostless "
                 "notion of its own: it marks the point source under the transient where SBW drops it as a point source, "
                 "which is the main split. This universe differs from baseline only in dp2_sample.")
    lines += ["## Numbers", ""]
    j = load_json("results/dp2_ddf_clean/dp2_cross_method_metrics.json")
    if not j:
        lines.append("- dp2_cross_method_metrics not readable here"); return lines
    a = j["subsets"]["all"]
    lines.append(f"- N {a['n']:,}; hostless: " + ", ".join(f"{k} {pct(v['frac'])}" for k, v in a["hostless"].items()))
    for pair, c in a["pairs"].items():
        both = c["same"]["k"] + c["different"]["k"]
        lines.append(f"- {pair}: same {c['same']['k']:,} ({pct(c['same']['frac'])}), different {c['different']['k']:,}, "
                     f"same given both hosted {pct(c['same']['k'] / both) if both else 'n/a'}; A-only hostless {c['a_only_hostless']['k']:,}, "
                     f"B-only hostless {c['b_only_hostless']['k']:,}, both hostless {c['both_hostless']['k']:,}")
    p = a["p_sbw_of_delight_host"]
    lines.append(f"- P_SBW of the DELIGHT host (N {p['n']:,}): > 0.5 {pct(p['gt_0p5']['frac'])}, 0.05..0.5 {pct(p['t_to_0p5']['frac'])}, < 0.05 {pct(p['lt_t']['frac'])}, "
                 f"SBW hostless {pct(p['sbw_hostless']['frac'])}; DELIGHT status of those: {p.get('delight_status')}")
    rm = j.get("reference_match_sensitivity", {})
    if rm:
        lines.append("- DELIGHT matched under each reference_match: " + ", ".join(f"{k} {pct(v['delight_matched']['frac'], 1)}" for k, v in rm.items()))
    lines.append(f"- median SN-to-DELIGHT separation {j.get('median_sep_sn_delight_arcsec', float('nan')):.2f} arcsec")
    return lines


def table2():
    lines = head("Table 2. Summary of every tier", [("baseline", "summary_table", "csv")],
                 "One row per tier on its headline denominator (rung 4 where a reference exists; the 690 rc_salt transients on DP2).",
                 "Columns: N detectable / host in the object table / host usable; mean P_SBW(reference) with standard error; "
                 "the four P bands; SBW argmax and DLR agreement with the reference (Wilson 68 %); DLR hostless; SBW-DLR "
                 "agreement; random positions rejected by each method; calibration error; the two fixed thresholds. DP2's "
                 "reference is the DLR host, so its agreement is same-host given DLR hosted and its calibration columns are empty.")
    lines += ["## Rows", ""]
    p = ROOT / "results/baseline/summary_table.csv"
    if readable(str(p.relative_to(ROOT))):
        rows = list(csv.DictReader(open(p)))
        keys = list(rows[0])
        lines.append("| " + " | ".join(keys) + " |")
        lines.append("|" + "---|" * len(keys))
        for r in rows:
            lines.append("| " + " | ".join(str(r[k])[:12] for k in keys) + " |")
    else:
        lines.append("- summary_table not readable here")
    return lines


def main():
    OUT.mkdir(exist_ok=True)
    for name, fn in (("fig1_placement", fig1), ("fig2_reliability", fig2), ("fig3_matrix_sims", fig3),
                     ("fig4_host_magnitude", fig4), ("fig5_hostless_roc", fig5), ("fig6_des_matrix", fig6),
                     ("fig7_dp2_matrix", fig7), ("fig8_dp2_consistency", fig8), ("fig9_delight", fig9),
                     ("table2_summary", table2)):
        try:
            text = "\n".join(fn()) + "\n"
        except Exception as e:   # a missing key in an older result file: say so instead of dying
            text = f"# {name}\n\ncaption generation failed: {type(e).__name__}: {e}\n"
            print(f"{name}: {e}", file=sys.stderr)
        (OUT / f"{name}.md").write_text(text)
        print(f"wrote captions/{name}.md")


if __name__ == "__main__":
    main()
