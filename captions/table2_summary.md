# Table 2. Summary of every tier

## Outputs
- `results/baseline/summary_table.csv` (content on this machine); made at commit `38033a0d` with `python src/summary_table.py --dp02 results/baseline/dp02_recovery.json --ou24 results/baseline/ou24_recovery.json --des results/baseline/des_recovery.json --dp2 results/baseline/dp2_disagreement_rate.json --output results/baseline/summary_table.csv`

## Sample

One row per tier on its headline denominator (rung 4 where a reference exists; the 690 rc_salt transients on DP2).

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

Columns: N detectable / host in the object table / host usable; mean P_SBW(reference) with standard error; the four P bands; SBW agreement with the reference (P_SBW(reference) >= 0.05, the reference host in the credible set) and DLR agreement (Wilson 68 %); DLR hostless; SBW-DLR agreement (the DLR host in SBW's credible set); random positions rejected by each method; calibration error; the two fixed thresholds. DP2's reference is the DLR host, so its SBW agreement is P_SBW(DLR host) >= 0.05 given DLR hosted and its calibration columns are empty.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Rows

| tier | reference | n_detectable | n_host_in_table | n_host_usable | mean_p_sbw_reference | mean_p_sbw_reference_se | p_gt_0p5 | p_t_to_0p5 | p_lt_t | sbw_hostless | sbw_agreement | sbw_agreement_ci68 | dlr_agreement | dlr_agreement_ci68 | dlr_hostless | sbw_dlr_agreement | random_rejected_sbw | random_rejected_dlr | n_random | calibration_error | calibration_n | p_true_threshold | sbw_hostless_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dp02 | DC2 truth host | 31,955 | 30,806 | 30,642 | 97.63 % | 0.0007 | 95.98 % | 1.04 % | 0.96 % | 2.02 % | 99.01 % | [98.95, 99.06] | 92.36 % | [92.21, 92.51] | 3.35 % | 95.48 % | 85.62 % | 85.18 % | 319,550 | 0.0227 | 30,023 | 0.0500 | 0.0100 |
| ou24 | OU2024 truth host | 19,058 | 17,033 | 15,851 | 98.64 % | 0.0008 | 95.92 % | 0.21 % | 0.15 % | 3.72 % | 99.05 % | [98.97, 99.12] | 97.04 % | [96.90, 97.17] | 2.18 % | 96.74 % | 91.91 % | 91.62 % | 190,580 | 0.0032 | 15,261 | 0.0500 | 0.0100 |
| des | DES-SN5YR published host | 1,623 | 1,609 | 1,594 | 99.04 % | 0.0014 | 99.31 % | 0.38 % | 0.00 % | 0.31 % | 99.94 % | [99.84, 99.98] | 99.37 % | [99.14, 99.54] | 0.44 % | 99.56 % | 75.77 % | 87.90 % | 16,230 | 0.0124 | 1,589 | 0.0500 | 0.0100 |
| dp2 | DLR host (no truth) | 690 |  |  | 98.70 % |  | 83.91 % | 0.14 % | 0.29 % | 14.06 % | 99.15 % |  |  |  | 14.93 % | 84.35 % | 83.80 % | 83.75 % | 6,900 |  |  |  |  |
