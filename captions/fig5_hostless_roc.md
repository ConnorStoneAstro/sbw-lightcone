# Fig 5. Hostless indicators on every tier (ROC)

## Outputs
- `results/baseline/hostless_roc_tiers.pdf` (annex pointer here; content on the machine that made it); made at commit `?` with `?`
- `results/baseline/dp02_hostless_metrics.json` (annex pointer here; content on the machine that made it); made at commit `6cfae476` with `python src/dp02_fig.py hostless_metrics --tier dp02 --associations results/baseline/dp02_associations.parquet --random results/baseline/dp02_random_associations.parquet --output results/baseline/dp02_hostless_metrics.json`
- `results/baseline/ou24_hostless_metrics.json` (annex pointer here; content on the machine that made it); made at commit `6cfae476` with `python src/dp02_fig.py hostless_metrics --tier ou24 --associations results/baseline/ou24_associations.parquet --random results/baseline/ou24_random_associations.parquet --output results/baseline/ou24_hostless_metrics.json`
- `results/baseline/des_hostless_metrics.json` (content on this machine); made at commit `17d075e3` with `python src/des_run.py fig hostless_metrics --tier des --associations results/baseline/des_associations.parquet --random results/baseline/des_random_associations.parquet --output results/baseline/des_hostless_metrics.json`
- `results/baseline/dp2_hostless_metrics.json` (annex pointer here; content on the machine that made it); made at commit `24c6577a` with `python src/dp02_fig.py hostless_metrics --tier dp2 --rung none --associations results/baseline/dp2_associations.parquet --random results/baseline/dp2_random_associations.parquet --output results/baseline/dp2_hostless_metrics.json`

## Sample

Real transients on rung 4 where a reference exists (DP0.2, OU2024, DES), all 690 rc_salt transients on DP2; random positions are ten per transient in the 45-90 arcsec annulus around it on every tier (2 Nov 2026; before, the designs differed by tier and so did the random baselines).

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

For SBW the hostless statistic is frac_light_beyond of the most probable host (flb_best): the fraction of that host's light that lies farther from its centre than the transient; for DLR it is the smallest d_DLR. Each curve sweeps the threshold; filled circles are the fixed operating points (flb 0.01, d_DLR 4) and open diamonds the Youden-J optima of the DP0.2 tier of the same universe (read from dp02_hostless_metrics, not frozen constants), so the other three tiers are held-out tests of whether a threshold tuned on DC2 transfers. Random positions get a P > 0.5 most probable host on every tier almost always, so hostlessness must come from flb, never from P.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- DC2 / DP0.2: result file not readable here
- OpenUniverse2024: result file not readable here
- DES-SN5YR: N real 1,594 (sample 1,623), random 16,230; at flb 0.01 kept 99.7 %, random rejected 75.8 %, J 0.755; at d_DLR 4 kept 99.6 %, rejected 87.9 %, J 0.875; own J optimum flb 0.133 (J 0.869), d_DLR 2.6 (J 0.910); DP0.2 optima used: flb 0.0215, d_DLR 2.7 (source ?); random positions landing on a P > 0.5 host 98.0 %
- Rubin DP2 rc_salt: result file not readable here
