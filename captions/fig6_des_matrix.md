# Fig 6. DES-SN5YR: SBW and DLR against the published host

## Outputs
- `results/baseline/des_agreement_matrix.pdf` (content on this machine); made at commit `17d075e3` with `python src/des_run.py fig agreement_matrix --tier des --associations results/baseline/des_associations.parquet --random results/baseline/des_random_associations.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/des_agreement_matrix.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005
- `results/baseline/des_recovery.json` (content on this machine); made at commit `17d075e3` with `python src/des_run.py recovery --associations results/baseline/des_associations.parquet --random results/baseline/des_random_associations.parquet --sample results/baseline/des_sample.parquet --sweep results/baseline/des_association_sweep.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/des_recovery.json`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

The 1,623 DES-SN5YR SNe Ia with a published host; the matrix is on rung 4 (published host matched to a DES DR2 object within 1 arcsec, kept by the object-table cuts, with a DLR ellipse). The reference is itself a DLR assignment on deeper SN-field coadds, so this is agreement, not recovery.

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

Same layout as Fig 3 with the published host as reference; candidates are DES DR2 objects within 30 arcsec with SExtractor light models; random positions ten per SN in the annulus. HOSTGAL2 (the published second host) explains most of the remaining 'other object' picks; the SNe whose published host is absent from DR2 or classed as a point source are the catalogue part of the ladder and appear only in the numbers.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- settings: profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness = ext05, candidate_quality = none
- ladder: {'rung2': 1623, 'rung3': 1609, 'rung4': 1594}; rung 4 N = 1,594
- SBW most probable = published 97.78 % [97.39, 98.12]; hosted and = published 97.54 %; DLR = published 97.60 % [97.19, 97.95]; mean P_SBW(published) 99.04 %
- published host status over the sample: {'absent': 14, 'point_source': 15, 'unusable': 0, 'quality': 0, 'kept': 1594}; matched in DR2 99.14 %
- HOSTGAL2: present 710, matched 588; of SBW picks that differ from the published host, HOSTGAL2 44.44 %; DLR 71.43 %
- matrix cells: p_gt_0p5/dlr_true 1,582, p_gt_0p5/dlr_wrong 0, p_gt_0p5/dlr_hostless 1, p_t_to_0p5/dlr_true 2, p_t_to_0p5/dlr_wrong 3, p_t_to_0p5/dlr_hostless 1, p_lt_t/dlr_true 0, p_lt_t/dlr_wrong 0, p_lt_t/dlr_hostless 0, sbw_hostless/dlr_true 0, sbw_hostless/dlr_wrong 0, sbw_hostless/dlr_hostless 5; random N 16,230, SBW rejects 75.77 %, DLR 87.90 %
- calibration error 0.0124 on 1,589 SNe
- per band/profile (same denominator as the headline): r/exponential SBW argmax 97.78 %, DLR 97.84 %
