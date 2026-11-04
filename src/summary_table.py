"""Table 2 of the paper: one row per tier with the same quantities under the same names, read from the
recovery / disagreement metric files (rung-4 headline where a reference exists; DP2 against the DLR host)."""
from __future__ import annotations

import argparse
import csv
import json


def pick(d, *keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return default
    return d


def truth_row(tier, reference, r):
    bands = r.get("p_true_bands", {})
    rm = r.get("random_matrix", {})
    return {"tier": tier, "reference": reference,
            "n_detectable": pick(r, "ladder", "rung2"), "n_host_in_table": pick(r, "ladder", "rung3"),
            "n_host_usable": pick(r, "ladder", "rung4"),
            "mean_p_sbw_reference": r.get("mean_p_sbw_true"), "mean_p_sbw_reference_se": r.get("mean_p_sbw_true_se"),
            "p_gt_0p5": bands.get("p_gt_0p5"), "p_t_to_0p5": bands.get("p_t_to_0p5"), "p_lt_t": bands.get("p_lt_t"),
            "sbw_hostless": bands.get("sbw_hostless"),
            "sbw_agreement": r.get("sbw_recovery"), "sbw_agreement_ci68": r.get("sbw_recovery_ci68"),
            "dlr_agreement": r.get("dlr_recovery"), "dlr_agreement_ci68": r.get("dlr_recovery_ci68"),
            "dlr_hostless": r.get("dlr_hostless_fraction"), "sbw_dlr_agreement": r.get("agreement"),
            "random_rejected_sbw": rm.get("sbw_rejected"), "random_rejected_dlr": rm.get("dlr_rejected"),
            "n_random": rm.get("n"), "calibration_error": r.get("expected_calibration_error"),
            "calibration_n": r.get("calibration_n_sne"), "p_true_threshold": r.get("p_true_threshold"),
            "sbw_hostless_cut": r.get("sbw_hostless_cut")}


def dp2_row(r):
    h, rnd = r["headline"], r["random_positions"]
    m = r.get("p_band_matrix", {})
    n = h["n"]
    hosted_dlr = 1 - h["dlr_hostless"]
    in_set = 1 - h["frac_p_sbw_of_dlr_host_lt_t"]   # DLR host in SBW's credible set, given DLR hosted
    return {"tier": "dp2", "reference": "DLR host (no truth)",
            "n_detectable": n, "n_host_in_table": None, "n_host_usable": None,
            "mean_p_sbw_reference": h.get("mean_p_sbw_of_dlr_host"), "mean_p_sbw_reference_se": None,
            "p_gt_0p5": m.get("p_dlr_gt_0p5/dlr_hosted", 0) / n if m else None,
            "p_t_to_0p5": m.get("p_dlr_t_to_0p5/dlr_hosted", 0) / n if m else None,
            "p_lt_t": m.get("p_dlr_lt_t/dlr_hosted", 0) / n if m else None,
            "sbw_hostless": h["sbw_hostless"],
            "sbw_agreement": in_set, "sbw_agreement_ci68": None,
            "dlr_agreement": None, "dlr_agreement_ci68": None,
            "dlr_hostless": h["dlr_hostless"], "sbw_dlr_agreement": in_set * hosted_dlr,
            "random_rejected_sbw": rnd["sbw_hostless"], "random_rejected_dlr": rnd["dlr_hostless"],
            "n_random": rnd["n"], "calibration_error": None, "calibration_n": None,
            "p_true_threshold": pick(r, "settings", "p_threshold"), "sbw_hostless_cut": pick(r, "settings", "sbw_hostless_cut")}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dp02", required=True); ap.add_argument("--ou24", required=True)
    ap.add_argument("--des", required=True); ap.add_argument("--dp2", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    rows = [truth_row("dp02", "DC2 truth host", json.load(open(a.dp02))),
            truth_row("ou24", "OU2024 truth host", json.load(open(a.ou24))),
            truth_row("des", "DES-SN5YR published host", json.load(open(a.des))),
            dp2_row(json.load(open(a.dp2)))]
    with open(a.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
