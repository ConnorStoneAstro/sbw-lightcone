# Fig 9. SBW, DLR and DELIGHT on the DP2 deep-drilling transients

## Outputs
- `results/dp2_ddf_clean/dp2_cross_method_matrix.pdf` (annex pointer here; content on the machine that made it); made at commit `?` with `?`
- `results/dp2_ddf_clean/dp2_cross_method_metrics.json` (annex pointer here; content on the machine that made it); made at commit `6fdffc57` with `python src/dp2_cross_method.py metrics --associations results/dp2_ddf_clean/dp2_associations.parquet --candidates results/dp2_ddf_clean/dp2_candidates.parquet --delight data/dp2_delight_yjkang_index.csv --reference-match pos_1arcsec --hostless-cut flb01 --p-threshold p005 --kang-ls-id-shift 20 --workers 16 --output results/dp2_ddf_clean/dp2_cross_method_metrics.json`
  - decisions in the recipe: reference_match = pos_1arcsec, sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

The cleaned DDF ParSNIP sample (dp2_sample = parsnip_ddf_clean: the DDF transients with a DELIGHT host from Y.-J. Kang, minus those on a pre-existing Legacy point source that is a Gaia star or a WISE AGN); DELIGHT's host is matched to a DP2 Object within 1 arcsec (reference_match = pos_1arcsec).

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

For each pair of methods a 2x2 of hosted / hostless with the hosted-hosted cell split into same host / different host; the title gives same host given both hosted. The last panel is the P_SBW band of the DELIGHT host over the transients whose DELIGHT position matched a candidate. DELIGHT has no hostless notion of its own: it marks the point source under the transient where SBW drops it as a point source, which is the main split. This universe differs from baseline only in dp2_sample.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- dp2_cross_method_metrics not readable here
