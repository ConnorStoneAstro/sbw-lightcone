# Fig 8. Band and light-model consistency of the SBW host on real data (DP2)

## Outputs
- `results/baseline/dp2_band_consistency.pdf` (annex pointer here; content on the machine that made it); made at commit `38033a0d` with `python src/dp2_fig.py band_consistency --associations results/baseline/dp2_associations.parquet --sweep results/baseline/dp2_association_sweep.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/dp2_band_consistency.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005
- `results/baseline/dp2_disagreement_rate.json` (content on this machine); made at commit `38033a0d` with `python src/dp2_fig.py disagreement_rate --associations results/baseline/dp2_associations.parquet --sweep results/baseline/dp2_association_sweep.parquet --random results/baseline/dp2_random_associations.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/dp2_disagreement_rate.json`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

The 690 rc_salt DP2 transients, the association recomputed in every band u-y and for every light model (exponential fit, free-n Sersic fit, moment Gaussian) from dp2_association_sweep.

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

(a) per band, how often the SBW most probable host and the DLR host differ from the r-band ones when both name one; (b) the SBW-DLR categories per band; (c) hostless fractions per band and light model; (d) SBW host agreement between light models per band. Without truth this is the stability of the answer, not its correctness; u is the exception for lack of flux.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- SBW host changes vs r band (exponential), given hosted in both: u 1.69 %, g 0.34 %, r 0.00 %, i 0.34 %, z 0.34 %, y 0.34 %
- hostless status changes vs r band: u 26.09 %, g 1.01 %, r 0.00 %, i 0.00 %, z 0.29 %, y 0.72 %
- same SBW host between light models (r band, both hosted): exponential_vs_sersic 100.00 % (N 591), exponential_vs_moments 99.82 % (N 556), sersic_vs_moments 99.82 % (N 557)
