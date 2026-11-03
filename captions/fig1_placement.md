# Fig 1. How SNe are placed with a host (DP0.2 and OU2024)

## Outputs
- `results/baseline/dp02_miss_accounting_fig.pdf` (content on this machine); made at commit `17c8eb80` with `python src/dp02_miss.py figure --legend outside --accounting results/baseline/dp02_miss_accounting.json --output results/baseline/dp02_miss_accounting_fig.pdf`
- `results/baseline/ou24_miss_accounting_fig.pdf` (content on this machine); made at commit `17c8eb80` with `python src/dp02_miss.py figure --legend outside --accounting results/baseline/ou24_miss_accounting.json --output results/baseline/ou24_miss_accounting_fig.pdf`
- `results/baseline/dp02_miss_accounting.json` (content on this machine); made at commit `b5df709b` with `python src/dp02_miss.py compute --sample results/baseline/dp02_sample.parquet --associations results/baseline/dp02_associations.parquet --candidates results/baseline/dp02_candidates.parquet --profile exponential --band r --search-radius r30 --dlr-max d4 --dlr-radius moments --extendedness ext05 --hostless-cut flb01 --p-threshold p005 --tier dp02 --workers 16 --output results/baseline/dp02_miss_accounting.json`
  - decisions in the recipe: sbw_profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness_cut = ext05, sbw_hostless_cut = flb01, p_true_threshold = p005
- `results/baseline/ou24_miss_accounting.json` (content on this machine); made at commit `b5df709b` with `python src/dp02_miss.py compute --sample results/baseline/ou24_sample.parquet --associations results/baseline/ou24_associations.parquet --candidates results/baseline/ou24_candidates.parquet --profile exponential --band r --search-radius r30 --dlr-max d4 --dlr-radius moments --extendedness ext05 --hostless-cut flb01 --p-threshold p005 --tier ou24 --workers 16 --output results/baseline/ou24_miss_accounting.json`
  - decisions in the recipe: sbw_profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness_cut = ext05, sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

Every detectable SN of the simulation (the sample); this is the only figure that shows the whole sample, so the reader sees the catalogue's losses and the method's losses side by side.

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

One stacked bar per method. From the bottom: hosted with the true host (green); hosted with a wrong host although the true host was available (red: the method's real misses; for SBW the most probable host is wrong, and the json counts how many of those still hold the true host at P >= 0.05); hosted with a wrong host because the true host was not in the object table (dark red: a catalogue effect, wrong by construction); rejected although the true host was available (light blue: the hostless cut fired); rejected with the true host not in the object table (grey). Beside the bar each fraction of the sample and, for the two wrong-host outcomes, the fraction of the *hosted* SNe, which is the contamination a user of the hosted sample meets. Under the bar: hosted fraction and contamination. The transients dropped as undetectable (peak r >= 24 or outside the imaged years) are not in the figure: their count is the survey selection function and is quoted from the ladder below. The choice to keep host-unmatched SNe in the sample (31 Oct 2026) is what makes the dark-red and grey segments possible: before, those SNe were dropped before association and nobody knew what the pipeline did with them.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

### DC2 / DP0.2 (settings: profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness = ext05, candidate_quality = none, truth_match = matches_truth, dp02_sample_cuts = physical, transient_peak_mag = r24)
- ladder: transients 394,022, hosted 394,021, with host row 354,053, detectable 39,352 (faint dropped 314,701), beyond the offset cut 7,397; sample 31,955 of which host-unmatched 1,149
- true host not in the object table: 1,313 (4.11 %): absent 1,149, removed by the cuts 164
- SBW: hosted 30,422 (95.20 %), rejected 1,533; true host 92.05 % of sample; wrong host, true host available 609 (1.91 % of sample, 2.00 % of hosted, of which 316 still hold the true host at P >= 0.05); wrong host, true host not in table 399 (1.31 % of hosted; absent 377, cut 22); rejected with true host available 619 (1.94 %); rejected, true host not in table 914; contamination 3.31 % [3.21, 3.42] of hosted
- DLR: hosted 29,994 (93.86 %), rejected 1,961; true host 88.57 % of sample; wrong host, true host available 1,314 (4.11 % of sample, 4.38 % of hosted); wrong host, true host not in table 378 (1.26 % of hosted; absent 367, cut 11); rejected with true host available 1,026 (3.21 %); rejected, true host not in table 935; contamination 5.64 % [5.51, 5.78] of hosted

### OpenUniverse2024 (settings: profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness = ext05, candidate_quality = none, truth_match = matches_truth, dp02_sample_cuts = physical, transient_peak_mag = r24)
- ladder: transients 511,108, hosted 511,108, with host row 510,928, detectable 21,383 (faint dropped 489,545), beyond the offset cut 2,325; sample 19,058 of which host-unmatched 2,025
- true host not in the object table: 3,207 (16.83 %): absent 2,025, removed by the cuts 1,182
- SBW: hosted 16,110 (84.53 %), rejected 2,948; true host 79.78 % of sample; wrong host, true host available 57 (0.30 % of sample, 0.35 % of hosted, of which 34 still hold the true host at P >= 0.05); wrong host, true host not in table 849 (5.27 % of hosted; absent 800, cut 49); rejected with true host available 590 (3.10 %); rejected, true host not in table 2,358; contamination 5.62 % [5.45, 5.81] of hosted
- DLR: hosted 16,341 (85.74 %), rejected 2,717; true host 80.71 % of sample; wrong host, true host available 123 (0.65 % of sample, 0.75 % of hosted); wrong host, true host not in table 836 (5.12 % of hosted; absent 795, cut 41); rejected with true host available 346 (1.82 %); rejected, true host not in table 2,371; contamination 5.87 % [5.69, 6.06] of hosted

