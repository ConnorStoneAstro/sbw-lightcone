# Fig 9. SBW, DLR and DELIGHT on the DP2 deep-drilling transients

## Outputs
- `results/baseline/dp2_cross_method_matrix.pdf` (content on this machine); made at commit `51552d94` with `python src/cross_method_matrix.py --metrics results/baseline/dp2_cross_method_metrics.json --subset all --output results/baseline/dp2_cross_method_matrix.pdf`
- `results/baseline/dp2_cross_method_metrics.json` (content on this machine); made at commit `51552d94` with `python src/dp2_cross_method.py metrics --associations results/baseline/dp2_ddf_associations.parquet --candidates results/baseline/dp2_ddf_candidates.parquet --delight data/dp2_delight_yjkang_index.csv --reference-match pos_1arcsec --hostless-cut flb01 --p-threshold p005 --kang-ls-id-shift 20 --workers 16 --output results/baseline/dp2_cross_method_metrics.json`
  - decisions in the recipe: reference_match = pos_1arcsec, sbw_hostless_cut = flb01, p_true_threshold = p005
- `results/baseline/dp2_ddf_sample.parquet` (annex pointer here; content on the machine that made it); made at commit `51552d94` with `python src/dp2_sample.py --index data/dp2_sn_candidates_index.csv --legacy-sweeps data/legacy_dr10_sweeps_index.csv --sample parsnip_ddf_clean --record-decisions --output results/baseline/dp2_ddf_sample.parquet`

## Sample

The cleaned DDF ParSNIP sample (dp2_ddf_sample, 2,601 transients: the DDF transients with a DELIGHT host from Y.-J. Kang, minus those on a pre-existing Legacy point source that is a Gaia star or a WISE AGN); DELIGHT's host is matched to a DP2 Object within 1 arcsec (reference_match = pos_1arcsec).

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

For each pair of methods a 2x2 of hosted / hostless with the hosted-hosted cell split into same host / different host; the title gives same host given both hosted. The last panel is the P_SBW band of the DELIGHT host over the transients whose DELIGHT position matched a candidate. DELIGHT has no hostless notion of its own: it marks the point source under the transient where SBW drops it as a point source, which is the main split. The DDF chain (dp2_ddf_sample, dp2_ddf_candidates, dp2_ddf_associations) is a fixed part of baseline (4 Nov 2026) made under the same decisions as the rc_salt chain; it serves this appendix figure only, since three quarters of the uncleaned DDF set sit on a star or AGN.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- N 2,601; hostless: SBW 23.26 %, DLR 27.99 %, DELIGHT 16.69 %, DELIGHT-SEx 3.31 %
- SBW|DLR: same 1,812 (69.67 %), different 24, same given both hosted 98.69 %; A-only hostless 37, B-only hostless 160, both hostless 568
- SBW|DELIGHT: same 1,419 (54.56 %), different 282, same given both hosted 83.42 %; A-only hostless 466, B-only hostless 295, both hostless 139
- SBW|DELIGHT-SEx: same 1,590 (61.13 %), different 350, same given both hosted 81.96 %; A-only hostless 575, B-only hostless 56, both hostless 30
- DLR|DELIGHT: same 1,411 (54.25 %), different 188, same given both hosted 88.24 %; A-only hostless 568, B-only hostless 274, both hostless 160
- DLR|DELIGHT-SEx: same 1,584 (60.90 %), different 239, same given both hosted 86.89 %; A-only hostless 692, B-only hostless 50, both hostless 36
- DELIGHT|DELIGHT-SEx: same 2,102 (80.82 %), different 42, same given both hosted 98.04 %; A-only hostless 371, B-only hostless 23, both hostless 63
- P_SBW of the DELIGHT host (N 2,167): > 0.5 65.48 %, 0.05..0.5 0.37 %, < 0.05 12.64 %, SBW hostless 21.50 %; DELIGHT status of those: {'kept': 1522, 'point_source': 645}
- DELIGHT matched under each reference_match: pos_1arcsec 83.3 %, pos_3arcsec 96.2 %, frac_semimajor 33.9 %
- median SN-to-DELIGHT separation 0.54 arcsec
