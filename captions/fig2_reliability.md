# Fig 2. Candidate-level reliability of P_SBW (DP0.2 and OU2024; DES in text)

## Outputs
- `results/baseline/dp02_calibration.pdf` (content on this machine); made at commit `6cfae476` with `python src/dp02_fig.py calibration --associations results/baseline/dp02_associations.parquet --hostless-cut flb01 --output results/baseline/dp02_calibration.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01
- `results/baseline/ou24_calibration.pdf` (content on this machine); made at commit `6cfae476` with `python src/dp02_fig.py calibration --tier ou24 --associations results/baseline/ou24_associations.parquet --hostless-cut flb01 --output results/baseline/ou24_calibration.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01
- `results/baseline/des_calibration.pdf` (content on this machine); made at commit `17d075e3` with `python src/des_run.py fig calibration --tier des --associations results/baseline/des_associations.parquet --hostless-cut flb01 --output results/baseline/des_calibration.pdf`
  - decisions in the recipe: sbw_hostless_cut = flb01

## Sample

Hosted SNe (SBW assigned a host at flb_best >= 0.01) whose true host is usable (rung 4), so some candidate is the right one; every candidate with P_SBW > 0.001 of those SNe enters, binned by its P_SBW. SNe that SBW rejects and SNe whose true host was not in the object table are left out: SBW's probabilities are not used on the first and cannot be right on the second.

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

Per bin, the fraction of candidates that are the true host against the bin's mean P_SBW; a calibrated probability sits on the 1:1 line in every bin, not only for the most probable host. Log-log over the full range (left), linear above 0.5 (right), candidate and true-host counts per bin below. The calibration error is the probability mass per SN that lands in the wrong bin: sum over bins of |true hosts - summed P| / N. The earlier top-host-only curve conditioned on the ranking and hid the low-P bins where SBW decides how much mass faint neighbours get (30 Oct 2026). The true-host probability is rounded to float32 like the candidate list so both histograms share the same bin membership.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- DC2 / DP0.2: calibration error 0.0227 per SN on 30,023 SNe; rung 4 N = 30,642; mean P_SBW(reference) 97.63 % +- 0.07; P bands: > 0.5 95.98 %, 0.05..0.5 1.04 %, < 0.05 0.96 %, SBW hostless 2.02 %
- OpenUniverse2024: calibration error 0.0032 per SN on 15,261 SNe; rung 4 N = 15,851; mean P_SBW(reference) 98.64 % +- 0.08; P bands: > 0.5 95.92 %, 0.05..0.5 0.21 %, < 0.05 0.15 %, SBW hostless 3.72 %
- DES-SN5YR (reference = published host): calibration error 0.0124 per SN on 1,589 SNe; rung 4 N = 1,594; mean P_SBW(reference) 99.04 % +- 0.14; P bands: > 0.5 99.31 %, 0.05..0.5 0.38 %, < 0.05 0.00 %, SBW hostless 0.31 %
