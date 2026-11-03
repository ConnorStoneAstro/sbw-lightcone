# Fig 1. How SNe are placed with a host (DP0.2 and OU2024)

## Outputs
- `results/baseline/dp02_miss_accounting_fig.pdf` (content on this machine); made at commit `6cfae476` with `python src/dp02_miss.py figure --accounting results/baseline/dp02_miss_accounting.json --output results/baseline/dp02_miss_accounting_fig.pdf`
- `results/baseline/ou24_miss_accounting_fig.pdf` (content on this machine); made at commit `6cfae476` with `python src/dp02_miss.py figure --accounting results/baseline/ou24_miss_accounting.json --output results/baseline/ou24_miss_accounting_fig.pdf`
- `results/baseline/dp02_miss_accounting.json` (content on this machine); made at commit `6cfae476` with `python src/dp02_miss.py compute --sample results/baseline/dp02_sample.parquet --associations results/baseline/dp02_associations.parquet --candidates results/baseline/dp02_candidates.parquet --profile exponential --band r --search-radius r30 --dlr-max d4 --dlr-radius moments --extendedness ext05 --hostless-cut flb01 --p-threshold p005 --tier dp02 --workers 16 --output results/baseline/dp02_miss_accounting.json`
  - decisions in the recipe: sbw_profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness_cut = ext05, sbw_hostless_cut = flb01, p_true_threshold = p005
- `results/baseline/ou24_miss_accounting.json` (content on this machine); made at commit `6cfae476` with `python src/dp02_miss.py compute --sample results/baseline/ou24_sample.parquet --associations results/baseline/ou24_associations.parquet --candidates results/baseline/ou24_candidates.parquet --profile exponential --band r --search-radius r30 --dlr-max d4 --dlr-radius moments --extendedness ext05 --hostless-cut flb01 --p-threshold p005 --tier ou24 --workers 16 --output results/baseline/ou24_miss_accounting.json`
  - decisions in the recipe: sbw_profile = exponential, band = r, search_radius = r30, dlr_max = d4, dlr_radius = moments, extendedness_cut = ext05, sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

Every detectable SN of the simulation (the sample); this is the only figure that shows the whole sample, so the reader sees the catalogue's losses and the method's losses side by side.

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

One stacked bar per method. From the bottom: hosted with the true host (green); hosted with a wrong host although the true host was available (red: the method's real misses; for SBW the most probable host is wrong, and the json counts how many of those still hold the true host at P >= 0.05); hosted with a wrong host because the true host was not in the object table (dark red: a catalogue effect, wrong by construction); rejected although the true host was available (light blue: the hostless cut fired); rejected with the true host not in the object table (grey). Beside the bar each fraction of the sample and, for the two wrong-host outcomes, the fraction of the *hosted* SNe, which is the contamination a user of the hosted sample meets. Under the bar: hosted fraction and contamination. The transients dropped as undetectable (peak r >= 24 or outside the imaged years) are not in the figure: their count is the survey selection function and is quoted from the ladder below. The choice to keep host-unmatched SNe in the sample (31 Oct 2026) is what makes the dark-red and grey segments possible: before, those SNe were dropped before association and nobody knew what the pipeline did with them.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- DC2 / DP0.2: result file not readable here
- OpenUniverse2024: result file not readable here
