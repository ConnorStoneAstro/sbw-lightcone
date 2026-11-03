# Fig 7. Rubin DP2 without a reference: SBW against DLR

## Outputs
- `results/baseline/dp2_agreement_matrix.pdf` (annex pointer here; content on the machine that made it); made at commit `6cfae476` with `python src/dp2_fig.py agreement_matrix --associations results/baseline/dp2_associations.parquet --random results/baseline/dp2_random_associations.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/dp2_agreement_matrix.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005
- `results/baseline/dp2_disagreement_rate.json` (content on this machine); made at commit `6cfae476` with `python src/dp2_fig.py disagreement_rate --associations results/baseline/dp2_associations.parquet --sweep results/baseline/dp2_association_sweep.parquet --random results/baseline/dp2_random_associations.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/dp2_disagreement_rate.json`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

The 690 rc_salt transients (SALT3-fit SNe Ia from the DP2 release candidate) with DP2 Object candidates within 30 arcsec and MultiProFit light models; random positions ten per transient in the annulus.

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

Rows: P_SBW of the DLR host > 0.5 / in [0.05, 0.5] / < 0.05, or SBW hostless; columns: the DLR host is SBW's most probable host / another object / DLR hostless. Without truth the question is whether the two methods name the same object when both assign one, and where their hostless calls differ (compact faint hosts whose moment ellipse is tiny).

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- settings: profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness = ext05, candidate_quality = none, exclude_suspect = keep
- N 690: same host 83.91 %, different host 0.43 %, SBW hostless only 0.72 %, DLR hostless only 1.59 %, both hostless 13.33 %; SBW hostless 14.06 %, DLR hostless 14.93 %; same host given DLR hosted 98.64 %; mean P_SBW(DLR host) 98.70 %; P(DLR host) < 0.05 0.85 %
- random positions N 6,900: SBW hostless 83.04 %, DLR hostless 83.16 %
