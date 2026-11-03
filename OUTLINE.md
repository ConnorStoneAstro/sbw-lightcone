# SBW host association: paper outline

Status, 2 Nov 2026. This file is the source the paper is written from:
the one configuration everything stands on, the vocabulary, the sections
with their figures (output ids, result paths, caption files), the
headline numbers, and what is still open. `index.md` stays the MyST
report. Every number here is read from a result file; the long-form
caption beside each figure (`captions/*.md`, written by
`python src/write_captions.py` from the result files after a
materialization) carries the full set of numbers and the decisions behind
that figure, so the text can be written from those files alone.

Framing (fixed): SBW is the method. The paper calibrates and characterises
it against truth (DC2, OU2024) and against references (DES-SN5YR published
hosts, DELIGHT on DP2), and relates it to DLR. Simulations give
*recovery*; real data gives *agreement*. "Which method wins" is never the
question; "is P_SBW honest, where does a host go at each step, how
contaminated is the hosted sample, and does a threshold tuned on
simulation transfer" are.

---

## 1. The configuration (settled)

The paper universe is `baseline`, which is exactly the spec defaults
(`universes/baseline.yaml`; since 30 Oct 2026). Every other universe is
baseline plus one change and appears only as a stated robustness check.

| decision | paper value | one-line reason |
|---|---|---|
| `sbw_profile` | exponential | the three light models are an output (profile sweep), not a universe |
| `band` | r | band dependence is an output |
| `search_radius` | r30 | d_DLR <= 4 already bounds DLR; r15 (DES-faithful) dropped from the story |
| `dlr_max` | d4 | the literature's nominal value; d5 never materialized |
| `dlr_radius` | moments | every published DLR uses a moment or Kron ellipse; with `half_light` the DLR ellipse is SBW's own, so the comparison is only informative with an independent DLR (appendix universe `physical` = half-light) |
| `extendedness_cut` | ext05 | the single largest catalogue effect in three tiers; shown as "true host not in the object table" everywhere and once as `ext_none` (OU2024) |
| `sbw_hostless_cut` | flb01 | fixed on DP0.2, never retuned; the ROC figure carries every tier's own optimum |
| `candidate_quality` | none | `basic` changes nothing, `strict` removes 10 % of true hosts; one sentence |
| `p_true_threshold` | p005 | the credible-set knee sits at 0.02 to 0.05 on both simulations |
| `exclude_suspect` | keep | never varied; not a paper decision |
| `truth_match` | matches_truth | `truth_nearest` differs for 75 of 306k; one sentence |
| `transient_peak_mag` | r24 | (new 31 Oct 2026) only transients an analyst would have: true peak r < 24, DC2 also peaking inside the five imaged years; the dropped fraction is the selection function, a number from the sample ladder, never a bar |
| `dp02_sample_cuts` | physical | SNe beyond 3 half-light radii of the host along its truth ellipse dropped; the only cut with a physical justification (DC2's untruncated n = 4 wings); host-unmatched SNe stay in the sample (31 Oct 2026) |
| `dp2_sample` | rc_salt (headline); parsnip_ddf_clean for DELIGHT only | rc_salt is the credible-SN sample; the DDF ParSNIP set is the only one DELIGHT covers and must be cleaned of stars and AGN first |
| `dp2_candidate_catalogue` | dp2_object | Legacy candidates are one robustness paragraph (hosted fraction drops ~15 points from depth) |
| `reference_match` | pos_1arcsec | DES 99.1 % matched; DELIGHT 83 % at 1", 96 % at 3" (one sentence) |

Removed: the DELIGHT training-set tier (the Förster set is not public);
`delight_*` outputs, `delight_catalogue`, `r60`. The paper has four
tiers: DC2 / DP0.2, OpenUniverse2024, DES-SN5YR, Rubin DP2.

Conventions fixed everywhere and never retuned: SBW hostless when
flb_best < 0.01 (or no SBW host); the true host counts as excluded when
P_SBW(true) < 0.05; DLR host = nearest in d_DLR with d_DLR <= 4 on the
second-moment ellipse, else hostless. Random positions: ten per transient,
uniform in area in the 45 to 90 arcsec annulus around it, on every tier
(2 Nov 2026).

---

## 2. Vocabulary: how a transient moves through the analysis

The reader who plans a cosmology analysis wants to know, for a transient
they have in hand, what the association does with it and how contaminated
the hosted sample is. The paper uses four words and one subset:

- **all**: every transient the tier provides (DC2 truth SNe; every OU2024
  transient; the 1,623 DES-SN5YR SNe Ia; the 690 DP2 rc_salt transients).
- **detectable** (the sample): the transients an analyst would have.
  Simulations: hosted in truth, with a host truth row, true peak r < 24
  (DC2 also peaking inside the five imaged years), within 3 half-light
  radii of the host along its truth ellipse. SNe whose true host has no
  catalogue object stay in it. Real tiers: every transient of the sample.
  The transients dropped here are the survey selection function: a number
  in the text (`sample_ladder` in the sample parquet metadata, copied into
  `*_miss_accounting.json["ladder"]`), never a bar.
- The pipeline then calls each detectable transient **rejected** (hostless
  call: SBW flb_best < 0.01 or no SBW host; DLR no candidate with d_DLR
  <= 4) or **hosted** (a host assigned).
- A hosted transient holds the **true host**, a **wrong host because the
  true host was not in the object table** (absent from the catalogue, or
  removed by the extendedness or usable-measurement cuts: a catalogue
  effect, wrong by construction), or a **wrong host although the true host
  was available** (a real miss; for SBW the json also counts how many of
  those still hold the true host at P >= 0.05). A rejected transient
  splits the same way by whether the true host was available.
- **Contamination** = wrong hosts of either cause over the hosted
  transients. This is the number the reader wants.
- **True host usable** (rung 4 in the code): the detectable transients
  whose true host survived the object-table cuts and has a DLR ellipse.
  It is the denominator of every method figure (reliability, matrices,
  recovery-vs-X, ROC real curve, credible sets, offsets, band and profile
  sweeps), because only there can an association be right. Every caption
  names it and gives N. `recovery.json` carries the counts at every rung
  (`ladder`) and the main numbers on each (`by_rung`).

Simulations give every row; real tiers skip the first split (they have no
truth selection) and DP2 has no reference at all, so its "hosted"
outcomes are agreement with DLR.

---

## 3. Paper outline with figures

Figure numbers are provisional. "(sim both)" means DP0.2 and OU2024 side
by side from the same script with `--tier`. Paths are in the `baseline`
universe unless stated.

### 1. Introduction
- Why host association matters for Rubin transients (redshifts, host
  properties, cosmology; Qu+24 delta w benchmark as the stakes).
- DLR as the standard and its known limits (moment ellipse, d_DLR <= 4
  nominal, hostless by construction, no probability).
- SBW in one paragraph: a probability over candidates from the host light
  at the transient position; `frac_light_beyond` as its natural hostless
  statistic. Point to the method paper / package.
- What this paper does: calibrate P_SBW against truth, account for what
  happens to every detectable transient, give the contamination of the
  hosted sample by cause, test whether thresholds fixed on one simulation
  transfer to another simulation and to two real surveys, and relate SBW
  to DLR and to DELIGHT. Not a horse race.

### 2. Method summary and conventions
- SBW: light model per candidate (exponential fit, r band), P over
  candidates within 30", `flb_best`. Same code on every tier; only the
  light-model adapter changes (CModel on DC2/OU2024, MultiProFit on DP2,
  SExtractor on DES DR2, Tractor on Legacy).
- DLR as implemented: moment ellipse, d_DLR <= 4, nearest in d_DLR.
- Object-table conventions that matter: extendedness >= 0.5 as the point
  source cut, the always-applied "no usable radius or flux" drop, DC2 bdE1
  bug so shapes come from PSF-subtracted HSM moments, Tractor and
  SExtractor angle conventions (checked on cutouts: `des_angle_check`,
  `dp2_legacy_angle_check`, text only).
- The vocabulary of Section 2 above, and recovery vs agreement; "SBW most
  probable host", never "SBW host = true".
- Thresholds fixed here and never retuned.
- Table 1: the four tiers (sample, candidate catalogue, light model,
  reference, N all / detectable / true host usable). Numbers from the
  ladders in `captions/fig1_placement.md`, `captions/fig6_des_matrix.md`,
  `captions/fig7_dp2_matrix.md`.

### 3. Data
- 3.1 DC2 / DP0.2: truth SNe, the SN truth summary for peak magnitudes,
  the DP0.2 visit window, MatchesTruth, the physical cut and why
  (untruncated Sersic wings; offset diagnostic in text), PA mirror as a
  one-sentence caveat with no consequence for association.
- 3.2 OpenUniverse2024: same products, w_2024_22, hand-built truth match,
  SNANA peak magnitudes, the host-mass bug (hosts ~1 dex too light), 15 s
  exposures. SB-weighted placement, so no wing tail.
- 3.3 DES-SN5YR: 1,623 Dovekie SNe Ia, DES DR2 candidates, published host
  is itself a DLR assignment on deeper coadds (so agreement, and the DLR
  column is expected to be near 100 %).
- 3.4 Rubin DP2: rc_salt 690 SALT3 candidates (headline); DDF ParSNIP
  7,823 for the DELIGHT comparison, cleaned to 2,601 by removing Legacy
  DR10 Gaia stars and WISE AGN under the transient (2,055 galaxy-only
  variant in text).
- 3.5 Random positions as the hostless control: ten per transient in the
  45 to 90 arcsec annulus, same cuts and scripts on every tier.

### 4. Calibration against truth (DP0.2 and OU2024)
- **Fig 1 (sim both): How SNe are placed with a host.**
  `results/baseline/dp02_miss_accounting_fig.pdf`,
  `results/baseline/ou24_miss_accounting_fig.pdf`; numbers in
  `results/baseline/{dp02,ou24}_miss_accounting.json` (`placement`,
  `ladder`); caption `captions/fig1_placement.md`. One stacked bar per
  method over the detectable SNe: hosted with the true host, the two
  wrong-host causes, the two rejected causes; the fraction of hosted SNe
  beside the wrong-host segments; hosted fraction and contamination under
  the bar. Message: SBW's real misses are a small fraction of a percent
  against DLR's few percent; the contamination of a hosted sample is set
  by the catalogue, not the method, and the hostless cut is what keeps
  it small; the wrong-host fraction inside "true host not in the object
  table" is nearly the same for both methods because, with the true host
  gone, the outcome depends on whether any other object lies in the
  search radius.
- **Fig 2 (sim both): candidate-level reliability of P_SBW.**
  `results/baseline/{dp02,ou24}_calibration.pdf` (DES in text or as a
  third panel: `results/baseline/des_calibration.pdf`); numbers
  `expected_calibration_error`, `calibration_n_sne` and `p_true_bands`
  in `*_recovery.json`; caption `captions/fig2_reliability.md`. Sample:
  hosted SNe whose true host is usable, every candidate with P > 0.001.
  Message: OU2024 is on the line in every bin; DC2 is under-confident
  below P ~ 0.1 (candidates given that little mass are the true host
  several times more often than stated), a feature of DC2's wing
  placement that the old top-host curve hid. Text: band and profile
  errors (`*_calibration_by_band`, `*_profile_sweep`).
- **Fig 3 (sim both): P_SBW(true) against DLR on SNe whose true host is
  usable.** `results/baseline/{dp02,ou24}_agreement_matrix.pdf`; cells in
  `*_recovery.json["agreement_matrix"]` and `["random_matrix"]`; caption
  `captions/fig3_matrix_sims.md`. Message: when DLR is wrong SBW usually
  still holds the true host at P > 0.5; random-position 2x2 beside it.
- **Fig 4 (sim both): recovery at fixed host magnitude.**
  `results/baseline/ou24_dp02_comparison.pdf`; numbers in
  `results/baseline/ou24_recovery_by_host_mag.json` (`rung4`, `rung2`);
  caption `captions/fig4_host_magnitude.md`. Message: at fixed host
  magnitude the two simulations agree; OU2024's lower headline is its
  faint hosts and the extendedness cut; the grey line is the catalogue's
  selection function per host magnitude.
- Text only: credible sets (knee and coverage from
  `*_credible_set_metrics.json`), profile sweep, band dependence, PA
  mirror, truth-match robustness, candidate-quality cuts (no help),
  the offset diagnostic.

### 5. Hostless calibration and the transfer of thresholds
- **Fig 5: hostless ROC on every tier.**
  `results/baseline/hostless_roc_tiers.pdf` (one panel per tier, drawn
  from the `curves` block of the four `*_hostless_metrics.json`); caption
  `captions/fig5_hostless_roc.md`. Message: SBW and DLR sit on similar
  curves; flb 0.01 and d_DLR 4 are the fixed points; the DP0.2 J optimum
  of the same universe (open diamonds) transfers to OU2024 at a small
  cost in J and moves on the real tiers because real SNe sit closer to
  their hosts; random positions still get a P > 0.5 most probable host
  almost always, so hostlessness must come from flb, never from P.
- Per-tier detail figures stay as outputs (`*_hostless_indicator.pdf`),
  not paper figures.

### 6. Agreement with references on real data
- 6.1 DES-SN5YR. **Fig 6: P-band matrix against the published host** on
  SNe whose published host is usable. `results/baseline/des_agreement_matrix.pdf`;
  numbers in `results/baseline/des_recovery.json`; caption
  `captions/fig6_des_matrix.md`. Message: SBW and DLR both agree with the
  published host at the 99 % level; the residuals are catalogue (absent
  or star-classed hosts) and HOSTGAL2 explains most "other" picks.
- 6.2 Rubin DP2 without a reference. **Fig 7: P-band matrix of the DLR
  host and hosted/hostless 2x2.** `results/baseline/dp2_agreement_matrix.pdf`;
  numbers in `results/baseline/dp2_disagreement_rate.json`; caption
  `captions/fig7_dp2_matrix.md`. Message: SBW and DLR name the same host
  whenever both assign one; the differences are hostless calls on compact
  faint hosts.
- **Fig 8: band and light-model consistency on real data.**
  `results/baseline/dp2_band_consistency.pdf`; numbers in
  `dp2_disagreement_rate.json` (`band_flip`, `profile_agreement`); caption
  `captions/fig8_dp2_consistency.md`.
- Text: `dp2_p_top_host` (numbers in `dp2_disagreement_rate.headline`),
  agreement vs d_DLR (collapses at d_DLR 4 by construction on every tier).

### 7. Comparison with DELIGHT on the DP2 deep-drilling fields
- Sample story in one paragraph: 7,823 ParSNIP transients, 74 % on a
  pre-existing Legacy point source, 79 % of those Gaia stars and 17 %
  WISE AGN; cleaned sample 2,601.
- **Fig 9: three-way SBW / DLR / DELIGHT.**
  `results/dp2_ddf_clean/dp2_cross_method_matrix.pdf`; numbers in
  `results/dp2_ddf_clean/dp2_cross_method_metrics.json`; caption
  `captions/fig9_delight.md`. Message: given both hosted, SBW and DELIGHT
  agree on the host for most transients; the split is a compact source
  under the transient (DELIGHT marks it, SBW drops it as a point source)
  versus a larger neighbour.
- Text: galaxy-only variant (`dp2_ddf_galaxy`), Legacy candidates
  (`dp2_ddf_clean_legacy`), P_SBW(DELIGHT object) bands, the
  reference-match sensitivity.

### 8. Discussion
- What P_SBW means operationally: the honest probability lets a user pick
  a credible set; one extra candidate buys about a point of coverage.
- Where a host goes: the catalogue (extendedness, deblending, depth) is
  where hosts are lost and where contamination comes from; the method's
  own misses are the small red segment of Fig 1. Implication for Rubin
  object tables: a lenient star cut or an extendedness-aware prior.
- Threshold transfer: flb 0.01 is a safe fixed point; DLR's tuned d_DLR
  does not transfer, 4 does.
- Limits: no stamp form yet; DC2 placement tail and its under-confident
  low-P bins; OU2024 host-mass bug; DES reference is DLR-derived; DELIGHT
  has no truth.
- Cosmology propagation if cheap (Qu+24 style) or explicitly deferred.

### 9. Summary
- One paragraph, and **Table 2** = `results/baseline/summary_table.csv`
  (caption `captions/table2_summary.md`): per tier N detectable / host in
  the object table / host usable, mean P_SBW(reference), the four P bands,
  argmax and DLR agreement, hostless fractions, random positions
  rejected, calibration error.

Appendix A: universes (`dlr_moments` = clean3 cut, `physical` =
half-light DLR, `ext_none`, `truth_nearest`, `sep10`, Legacy candidates,
galaxy-only DP2) as one table of headline numbers from their
`*_recovery.json`. Appendix B: the DC2 PA mirror and offset diagnostic.
Appendix C: Tractor and SExtractor conventions.

---

## 4. Headline numbers (baseline; filled from the result files)

See the `## Numbers` section of each caption file for the full set. This
table is the short version the text will quote; it is rewritten after
every materialization.

NUMBERS_TABLE

---

## 5. Inconsistencies found on 30 Oct 2026 and their status

- 4.1 Headline universe not the default: fixed 30 Oct (baseline = spec
  defaults; `physical` and `dlr_moments` keep their ids as the appendix
  universes with the half-light DLR and the clean3 cut respectively).
- 4.2 Two DLR definitions: fixed, moments everywhere; Legacy candidates
  have no moments and say so.
- 4.3 Denominators differ between figures: fixed 31 Oct (every method
  figure on "true host usable"; Fig 1 alone shows the whole sample) and
  2 Nov (DES per-band and agreement blocks on the same denominator as
  the headline).
- 4.4 Random-position designs differ by tier: fixed 2 Nov, ten per
  transient in the 45 to 90 arcsec annulus on every tier.
- 4.5 Two headline metrics: fixed by `summary_table` (mean P first, the
  P bands, argmax agreement) and the captions' vocabulary.
- 4.6 DP2 samples and universes: the paper shows rc_salt (Sections 5, 6)
  and dp2_ddf_clean (Section 7); the rest are sentences.
- 4.7 Thresholds tuned in the text: fixed operating points everywhere;
  the optima appear only on Fig 5 and in `*_hostless_metrics.json`.
- 4.8 DES r15: dropped; `des_faithful` is an appendix universe.
- 4.9 Extendedness on real data: one paragraph covering OU2024, DES and
  DP2 (the "true host not in the object table" outcome of Fig 1 is the
  same effect); `ext_none` on DP2 not materialized.
- 4.10 Stale spec text: DELIGHT tier removed; `summary_table` has a
  recipe; the miss figure is one panel.
- 4.11 DP0.2 optima hard-coded: fixed 2 Nov, every tier reads the DP0.2
  hostless and credible-set metrics of its own universe.
- 4.12 DES DLR 0.971 vs 0.973: fixed 2 Nov, the per-band table and the
  agreement block use the headline's denominator.
- 4.13 Settings blocks omit the sample decisions: fixed 2 Nov, samples
  and DP2 candidates record their decisions in parquet metadata and the
  association copies them into every settings block.
- 4.14 Tier vocabulary and CI formats: `summary_table` is the forcing
  function; captions quote Wilson 68 % intervals.
- 4.15 Smaller items: cat3 wording unified; `dp2_p_top_host` numbers in
  `dp2_disagreement_rate.headline`; grid statement in the ROC caption.
- 4.16 Stale numbers in AGENTS.md: the paper quotes only the caption
  files and result files; AGENTS.md keeps its history as history.

## 6. Still open

- `index.md` (the MyST report) has not been rewritten to the new
  framing.
- `ext_none` on DP2 and a DP2 `ext_none` sentence are not materialized.
- The `ou24_dp02_comparison` suptitle overlaps the panel titles
  (cosmetic).
- The paper-figure styling (fonts, sizes, colours) is the diagnostic
  styling; a final pass once the text exists.
