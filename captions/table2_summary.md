# Table 2. Summary of every tier

## Outputs
- `results/baseline/summary_table.csv` (annex pointer here; content on the machine that made it); made at commit `?` with `?`

## Sample

One row per tier on its headline denominator (rung 4 where a reference exists; the 690 rc_salt transients on DP2).

Vocabulary: *all* transients the tier provides; *detectable* = the sample (simulations: hosted, with a host truth row, true peak r < 24 and, for DC2, peaking inside the five imaged years, within 3 half-light radii of the host along its truth ellipse; real tiers: every transient of the sample); the pipeline then calls each detectable transient *rejected* (hostless) or *hosted*; a hosted transient holds the *true host*, a *wrong host because the true host was not in the object table* (absent, or removed by the extendedness / usable-measurement cuts), or a *wrong host although the true host was available* (a real miss). *True host usable* (rung 4) = the subset whose true host survived the object-table cuts and has a DLR ellipse; it is the denominator of every method figure.

## Decisions behind it

Columns: N detectable / host in the object table / host usable; mean P_SBW(reference) with standard error; the four P bands; SBW argmax and DLR agreement with the reference (Wilson 68 %); DLR hostless; SBW-DLR agreement; random positions rejected by each method; calibration error; the two fixed thresholds. DP2's reference is the DLR host, so its agreement is same-host given DLR hosted and its calibration columns are empty.

Fixed operating points, never retuned: SBW hostless when flb_best < 0.01 (or no SBW host); the true host counts as excluded when P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the second-moment ellipse, else DLR hostless.

## Rows

- summary_table not readable here
