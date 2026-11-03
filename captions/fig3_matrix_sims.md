# Fig 3. SBW probability of the true host against DLR, on SNe whose true host is usable (DP0.2 and OU2024)

## Outputs
- `results/baseline/dp02_agreement_matrix.pdf` (content on this machine); made at commit `6cfae476` with `python src/dp02_fig.py agreement_matrix --survivors-only --associations results/baseline/dp02_associations.parquet --random results/baseline/dp02_random_associations.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/dp02_agreement_matrix.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005
- `results/baseline/ou24_agreement_matrix.pdf` (content on this machine); made at commit `6cfae476` with `python src/dp02_fig.py agreement_matrix --survivors-only --tier ou24 --associations results/baseline/ou24_associations.parquet --random results/baseline/ou24_random_associations.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/ou24_agreement_matrix.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

Rung 4: SNe whose true host passed the object-table cuts for both methods, so every cell is an association outcome. The random-position panel beside it uses the ten-per-SN annulus draws under the same cuts.

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

Rows: P_SBW(true host) > 0.5, in [0.05, 0.5], < 0.05 (excluded), or SBW hostless at flb_best < 0.01; columns: DLR host = true, DLR another object, DLR hostless. Counts and fractions of rung 4 in every cell. The random 2x2 shows how often each method correctly rejects a position with no host. The 'both wrong, same object' tail of DC2 comes from its untruncated Sersic placement; the physical offset cut (3 R_e) removes most of it (dp02_offset_diagnostic).

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- DC2 / DP0.2: N = 30,642; p_gt_0p5/dlr_true 28,273 (92.27 %), p_gt_0p5/dlr_wrong 676 (2.21 %), p_gt_0p5/dlr_hostless 462 (1.51 %), p_t_to_0p5/dlr_true 17 (0.06 %), p_t_to_0p5/dlr_wrong 294 (0.96 %), p_t_to_0p5/dlr_hostless 8 (0.03 %), p_lt_t/dlr_true 2 (0.01 %), p_lt_t/dlr_wrong 286 (0.93 %), p_lt_t/dlr_hostless 5 (0.02 %), sbw_hostless/dlr_true 10 (0.03 %), sbw_hostless/dlr_wrong 58 (0.19 %), sbw_hostless/dlr_hostless 551 (1.80 %); SBW argmax = true 97.92 % [97.84, 98.00], DLR = true 92.36 % [92.21, 92.51], SBW most probable = DLR 94.12 %; random positions N = 386,248: SBW rejects 85.74 %, DLR rejects 85.07 %
- OpenUniverse2024: N = 15,851; p_gt_0p5/dlr_true 15,062 (95.02 %), p_gt_0p5/dlr_wrong 58 (0.37 %), p_gt_0p5/dlr_hostless 84 (0.53 %), p_t_to_0p5/dlr_true 2 (0.01 %), p_t_to_0p5/dlr_wrong 31 (0.20 %), p_t_to_0p5/dlr_hostless 1 (0.01 %), p_lt_t/dlr_true 6 (0.04 %), p_lt_t/dlr_wrong 17 (0.11 %), p_lt_t/dlr_hostless 0 (0.00 %), sbw_hostless/dlr_true 312 (1.97 %), sbw_hostless/dlr_wrong 17 (0.11 %), sbw_hostless/dlr_hostless 261 (1.65 %); SBW argmax = true 98.79 % [98.70, 98.87], DLR = true 97.04 % [96.90, 97.17], SBW most probable = DLR 96.46 %; random positions N = 521,766: SBW rejects 91.78 %, DLR rejects 91.58 %
