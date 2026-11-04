# Fig 4. Recovery at fixed host magnitude, DP0.2 against OU2024

## Outputs
- `results/baseline/ou24_dp02_comparison.pdf` (content on this machine); made at commit `b5df709b` with `python src/ou24_host_mag.py dp02_comparison --associations results/baseline/ou24_associations.parquet --dp02-associations results/baseline/dp02_associations.parquet --output results/baseline/ou24_dp02_comparison.pdf`
- `results/baseline/ou24_recovery_by_host_mag.json` (content on this machine); made at commit `38033a0d` with `python src/ou24_host_mag.py by_host_mag --associations results/baseline/ou24_associations.parquet --hostless-cut flb01 --p-threshold p005 --output results/baseline/ou24_recovery_by_host_mag.json`
  - decisions in the recipe: sbw_hostless_cut = flb01, p_true_threshold = p005

## Sample

Rung 4 of both simulations (true host usable), binned by the true host's r magnitude; a grey line per tier gives the fraction of sample hosts that reach rung 4 in each bin (the catalogue's selection function).

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

OU2024 hosts are drawn about 1 dex too light (OU2024 paper section 3.2.4), so the population-averaged numbers of the two simulations are not comparable; at fixed host magnitude the association is the same geometry and photometry problem in both. Curves: SBW mean P_SBW(true) and DLR fraction correct, DP0.2 solid, OU2024 dashed (all z, and z < 1 to match DC2's range); the third panel is the host-r distributions. The text box gives OU2024 reweighted to DP0.2's host-r distribution and both tiers at host r < 24.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Numbers

- OU2024 rung4 (true host usable by both methods), all_z: N 15,851, mean P 98.64 %, DLR 97.04 %; host r < 24: N 12,076, mean P 99.42 %, DLR 96.71 %
- OU2024 rung4 (true host usable by both methods), z_lt_1: N 15,766, mean P 98.65 %, DLR 97.04 %; host r < 24: N 12,064, mean P 99.42 %, DLR 96.71 %
- OU2024 rung2 (detectable (the sample)), all_z: N 19,058, mean P 82.04 %, DLR 80.71 %; host r < 24: N 12,643, mean P 94.97 %, DLR 92.38 %
- OU2024 rung2 (detectable (the sample)), z_lt_1: N 18,705, mean P 83.15 %, DLR 81.79 %; host r < 24: N 12,631, mean P 94.96 %, DLR 92.37 %
- DP0.2 rung 4: N 30,642, mean P 97.63 %, DLR 92.36 %
