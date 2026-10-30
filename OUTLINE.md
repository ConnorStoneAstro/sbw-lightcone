# SBW host association: paper outline and consistency audit

Status note, 30 Oct 2026. This is a planning document, not the report
(`index.md` stays the MyST report and will reference outputs by path). It
does three things: proposes the one set of choices the paper is written
under, lays out the sections with their figures, and lists where the
results made so far disagree with each other in sample, denominator,
threshold or definition. Numbers quoted here come from `AGENTS.md` and the
result JSONs as of main 8cc4bf2; every one of them is an argument for a
choice, not a number to copy into the text.

Framing (fixed): SBW is the method. The paper calibrates and characterises
it against truth (DC2, OU2024) and against references (DES-SN5YR published
hosts, DELIGHT on DP2), and relates it to DLR. Simulations give
*recovery*; real data gives *agreement*. "Which method wins" is never the
question; "is P_SBW honest, where does SBW lose the host, and does a
threshold tuned on simulation transfer" are.

---

## 1. One set of choices for the whole paper

Right now the headline numbers come from the `baseline` universe, which is
*not* the spec default (it keeps the 3 arcsec sample cut "for
comparability with the first run") and whose DLR is the package's
half-light ellipse rather than the literature's moment ellipse. Each tier
then has its own extra flavour (`des_faithful`, `dp2_ddf_clean`, …). The
reader should meet one configuration, with everything else as a stated
robustness check. Proposed paper configuration, decision by decision:

| decision | paper value | status today | verdict |
|---|---|---|---|
| `sbw_profile` | exponential | default, every universe | keep; the three forms are an *output* (profile sweep), not a universe |
| `band` | r | default | keep; band dependence is an output |
| `search_radius` | r30 | default; `des_faithful` uses r15 | keep r30 everywhere, drop r15 from the story (d_DLR <= 4 already bounds the DLR search; SBW does not need the DES radius) |
| `dlr_max` | d4 | default | keep; d5 never materialized, drop it or leave as an unmade universe |
| `dlr_radius` | **moments** | default is `half_light`; `dlr_moments` universe made for DP0.2, OU2024, DP2; DES only via `des_faithful` | **change the paper DLR to moments.** Every published DLR uses a moment or Kron ellipse; with `half_light` the DLR ellipse is the one SBW evaluates, so on DP2 the two methods agree "by construction" (AGENTS: 85.65 % same, 2 different). The SBW-vs-DLR comparison is only informative with an independent DLR. Cost: the Legacy-candidate universes have no moments (see 4.6). |
| `extendedness_cut` | ext05 | default | keep, but it is the single largest catalogue effect in three tiers (OU2024 true host lost 18 %; DP2 DELIGHT host dropped 67 %; DES published host classed as star 15/29 failures). Report it as cat2 everywhere and show `ext_none` once (OU2024) |
| `sbw_hostless_cut` | flb01 | default; J-optima 0.0215 (DP0.2), 0.010 (OU2024), 0.091 (DP2 rc_salt), 0.133 (DES) | keep 0.01 as the fixed operating point, chosen on DP0.2 and *not* retuned; the ROC figures carry the optima. Say explicitly that real SNe sit closer to their hosts than DC2's, so the optimum moves up on real data |
| `candidate_quality` | none | default | keep `none`; `basic` changes nothing (mean P 0.9705 to 0.9702) and `strict` removes 10 % of true hosts. Drop panel (c) of the miss-accounting figure; one sentence in the text |
| `p_true_threshold` | p005 | default; credible-set knee 0.034 (DP0.2), 0.046 (physical), 0.051 (OU2024) | keep 0.05; justified by the knee on both tiers |
| `exclude_suspect` | keep | default, never varied | keep; not a paper decision, consider removing from the spec |
| `truth_match` | matches_truth | default; `truth_nearest` checked (75 of 306k differ) | keep; one sentence that the same-wrong-object cell is not a labelling artefact |
| `dp02_sample_cuts` | **physical** | spec default is `physical`, but `baseline`, `dlr_moments`, `truth_nearest`, `ext_none` all run `clean3` | **make the paper sample `physical` on both simulated tiers.** It is the only cut with a physical justification (DC2's untruncated n = 4 wings; recovery flat to ~1.5 R_e then climbing); clean3 is a 3 arcsec cap that penalises large hosts and exists only for continuity with a run that no longer exists (pre-fallback-removal numbers are already superseded) |
| `dp2_sample` | rc_salt for SBW/DLR/hostless; parsnip_ddf_clean for DELIGHT | both made | keep two samples but say why once: rc_salt is the credible-SN sample (SALT3 fit), the DDF ParSNIP set is the only one DELIGHT covers and must be cleaned of stars/AGN first. Do not show hostless ROCs for both |
| `dp2_candidate_catalogue` | dp2_object | default | keep; Legacy candidates are a one-paragraph robustness check (hosted fraction drops ~15 points from depth), not a parallel set of figures |
| `reference_match` | pos_1arcsec | default | keep (DES 99.14 % matched; DELIGHT 90.4 %, 97.8 % at 3"); the 3" and frac_semimajor options need not appear |
| `delight_catalogue` | – | never materialized | the Förster training set is not public; **remove the DELIGHT tier** (`delight_*` outputs, `search_radius: r60`, the "five tiers" sentence in the analysis description) or mark it parked. The paper has four tiers |

Net effect: the paper universe is "spec defaults + `dlr_radius: moments`".
Done on 30 Oct 2026: `moments` is the spec default, `baseline.yaml` is the
spec defaults, and the appendix universes are `physical` (baseline with the
half-light DLR) and `dlr_moments` (baseline with clean3), kept under their
historical ids because a universe id is baked into every rendered recipe
and renaming them would make their current results stale. Nothing is
materialized yet; baseline and the one-change universes are stale by
decision change.

---

## 2. The sample ladder: one denominator convention per tier

Today the same tier reports three different N. DP0.2: 321,161 (hosted SNe,
miss-accounting denominator), 306,332 (sample after matcher and 3" cut,
hostless ROC and recovery), 304,454 (cat1+cat2 survivors, matrix). OU2024:
507,988 / 262,710 / 203,096. DES: 1,623 (matrix, ROC) vs 1,594 (survivors,
quoted in text). DP2: 690 vs 7,823 vs 2,601 vs 2,055. Each is defensible,
but the paper should walk one ladder per tier, print it once as a table,
and then say which rung every figure stands on.

Proposed ladder (same words in every tier):

1. **Transients**: every SN/transient the tier provides (DC2 truth SNe
   with z < 1 that are hosted; OU2024 all; DES 1,623; DP2 690 or 7,823).
2. **With a host row** (simulations only): the truth host exists and the
   SN sits within the physical cut. This is the *sample* and the
   denominator of the miss accounting (cat1 counted here).
3. **Host in catalogue** (cat1 passed): the reference host has a primary
   object within 1" (sims: truth match; DES: published host matched in
   DR2; DELIGHT: position matched to a DP2 object).
4. **Host usable** (cat2 passed): the host survives ext05 and has a radius
   and flux. This is the denominator of the agreement matrices, of
   calibration, and of recovery-vs-X curves.
5. Within rung 4, SBW either keeps the host (P > t), excludes it (P < t),
   or calls hostless (cat3/cat4).

Rule: every figure caption names its rung and N; the miss-accounting bar
chart is the only place rungs 2 to 5 appear together. Real-data tiers skip
rung 2 and say so.

The two real-data tiers currently stand on rung 1 for their matrices while
the simulations stand on rung 4. The DES matrix on 1,594 survivors (SBW
99.56 %, DLR 98.87 % faithful 99.37 %) is the like-for-like number; the
1,623 version then belongs in the miss accounting bar (cat1 14 unmatched,
cat2 15 stars).

---

## 3. Paper outline with figure placement

Figure numbers are provisional; "(sim both)" means one figure with DP0.2
and OU2024 side by side or stacked, made from the same script with
`--tier`, which is already how they are produced.

### 1. Introduction
- Why host association matters for Rubin transients (redshifts, host
  properties, cosmology; Qu+24 delta w benchmark as the stakes).
- DLR as the standard and its known limits (moment ellipse, d_DLR <= 4
  nominal, hostless by construction, no probability).
- SBW in one paragraph: a probability over candidates from the host light
  at the transient position; `frac_light_beyond` as its natural hostless
  statistic. Point to the method paper / package.
- What this paper does: calibrate P_SBW against truth, account for every
  lost host by cause, test whether thresholds fixed on one simulation
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
- The miss ladder (Section 2 above) and the vocabulary: recovery vs
  agreement; "SBW most probable host", never "SBW host = true".
- Thresholds fixed here and never retuned: flb < 0.01 hostless; P < 0.05
  excluded; d_DLR <= 4.
- Table 1: the four tiers (sample, candidate catalogue, light model,
  reference, N at each rung).

### 3. Data
- 3.1 DC2 / DP0.2: truth SNe, MatchesTruth, the physical cut and why
  (untruncated Sersic wings; offset diagnostic in text), PA mirror as a
  one-sentence caveat with no consequence for association.
- 3.2 OpenUniverse2024: same products, w_2024_22, hand-built truth match,
  the host-mass bug (hosts ~1 dex too light, median host r 25.1 vs 22.3),
  15 s exposures. SB-weighted placement, so no wing tail.
- 3.3 DES-SN5YR: 1,623 Dovekie SNe Ia, DES DR2 candidates, published host
  is itself a DLR assignment on deeper coadds (so agreement, and the
  DLR column is expected to be near 100 %).
- 3.4 Rubin DP2: rc_salt 690 SALT3 candidates (headline); DDF ParSNIP
  7,823 for the DELIGHT comparison, cleaned to 2,601 by removing Legacy
  DR10 Gaia stars and WISE AGN under the transient (and 2,055 galaxy-only
  variant in text).
- 3.5 Random positions as the hostless control: one draw per SN within
  the SN's own field unit; same cuts and same scripts (see 4.4 for the
  current mismatch).

### 4. Calibration against truth (DP0.2 and OU2024)
- **Fig 1 (sim both): miss accounting bars** (`dp02_miss_accounting_fig`,
  `ou24_miss_accounting_fig`, panels a+b only). Message: SBW loses hosts
  to the catalogue (cat1/cat2) and to its own hostless cut; genuine wrong
  picks are 0.3 to 0.7 %, against DLR's 3 to 5 %. On OU2024 half the
  hosts never reach the catalogue.
- **Fig 2 (sim both): reliability curve** (`dp02_calibration`,
  `ou24_calibration`). Message: P_SBW tracks the true-host fraction
  (ECE ~0.01 on DP0.2 physical); overconfident in the 0.5 to 0.95 band.
  Text: by band and by profile ECEs, not figures.
- **Fig 3 (sim both): agreement matrix on cat1+cat2 survivors**
  (`dp02_agreement_matrix`, `ou24_agreement_matrix`). Message: when DLR
  is wrong SBW usually still holds the true host at P > 0.5; the "both
  wrong, same object" cell is the sim's placement tail, not a method
  failure; random-position 2x2 beside it.
- **Fig 4 (sim both): recovery vs host magnitude**
  (`ou24_dp02_comparison`, which already overlays the two tiers and
  carries the host-r histograms). Message: at fixed host magnitude the
  two simulations agree to a point; OU2024's low headline is its faint
  hosts and the extendedness cut. Reweighted numbers in text.
- Text only: credible sets (knee 0.03 to 0.05, coverage 98.7 % at mean
  size 1.06 on DP0.2), profile sweep (de Vaucouleurs best on DP0.2,
  moments on OU2024, differences ~0.5 to 2 points), band dependence
  (u fails from missing shapes, g to y within 1 point on DP0.2), PA mirror,
  truth-match robustness, candidate-quality cuts (no help).

### 5. Hostless calibration and the transfer of thresholds
- **Fig 5 (four panels, one per tier): hostless ROC** with the DP0.2
  optimum marked (`*_hostless_indicator` panel d from each tier, or a new
  combined figure: `dp02`, `ou24`, `des`, `dp2` rc_salt). Message: SBW and
  DLR sit on the same ROC; flb 0.01 and d_DLR 4 are nearly the same
  ellipse; the flb cut chosen on DP0.2 costs < 0.01 in Youden J on OU2024
  and 0.02 to 0.04 on DP2; on DES the optimum moves to 0.13 because SNe
  sit at median d_DLR 0.67 from their hosts.
- Text: random positions still get a P > 0.5 host 99 % of the time on
  every tier, so hostlessness must come from flb, never from P.
- Possibly **Fig 5b**: panel (a) distribution of flb_best for real vs
  random on one tier, if the reader needs to see why.

### 6. Agreement with references on real data
- 6.1 DES-SN5YR. **Fig 6: P-band matrix vs published host** on survivors
  (`des_agreement_matrix` restricted as the sim matrices are). Message:
  SBW 99.6 %, DLR 99.4 % (moments); the 29 residuals are catalogue
  (unmatched or star), and HOSTGAL2 explains most of the remaining
  "other" picks. Reliability curve ECE 0.014 in text (or as a panel of
  Fig 2 if the figure takes three tiers).
- 6.2 Rubin DP2 without a reference. **Fig 7: P-band matrix of the DLR
  host and hosted/hostless 2x2** (`dp2_agreement_matrix`, rc_salt,
  moments DLR). Message: SBW and DLR name the same host whenever both
  assign one (98 %); the differences are hostless calls on compact faint
  hosts whose moment ellipse is tiny.
- **Fig 8: band and light-model consistency on real data**
  (`dp2_band_consistency`; DES version in text or as a second row).
  Message: the SBW host is stable across g to y (< 0.3 % flips) and across
  light models (100 % exponential vs Sersic when both hosted); u is the
  exception for lack of flux.
- Text: `dp2_p_top_host` (real-data SBW at least as decisive as on DC2),
  agreement vs d_DLR (collapses at d_DLR 4 by construction on every tier,
  which is why it need not be a figure).

### 7. Comparison with DELIGHT on the DP2 deep-drilling fields
- Sample story in one paragraph: 7,823 ParSNIP transients, 74 % on a
  pre-existing Legacy point source, 79 % of those Gaia stars and 17 %
  WISE AGN; cleaned sample 2,601.
- **Fig 9: three-way hosted/hostless and same-host matrix SBW / DLR /
  DELIGHT on the cleaned sample** (from `dp2_cross_method_metrics`,
  dp2_ddf_clean; needs a figure recipe, today it is a JSON). Message:
  given both hosted, SBW and DELIGHT agree 83 %; the remaining split is a
  compact source under the transient (DELIGHT marks it, SBW drops it as a
  point source) versus a larger neighbour.
- Text: galaxy-only variant (82 %), Legacy candidates (same agreement,
  15 points fewer hosted), P_SBW(DELIGHT object) > 0.5 for 65 %.
- Maybe a small gallery (4 to 6 cutouts) of the compact-source cases; the
  current galleries are diagnostic, not paper figures.

### 8. Discussion
- What P_SBW means operationally: the honest probability lets a user pick
  a credible set; one extra candidate buys ~1.2 points.
- Where SBW loses hosts: the catalogue (extendedness, deblending, depth),
  not the association. Implication for Rubin object tables: a lenient
  star cut or an extendedness-aware prior.
- Threshold transfer: flb 0.01 is a safe fixed point; DLR's tuned d_DLR
  2.7 does not transfer, 4 does.
- Limits: no stamp form yet; DC2 placement tail; OU2024 host-mass bug;
  DES reference is DLR-derived; DELIGHT has no truth.
- Cosmology propagation if cheap (Qu+24 style) or explicitly deferred.

### 9. Summary
- One paragraph; the summary table (`summary_table`, currently with no
  recipe) printed as Table 2: per tier, N at rung 4, mean P_SBW(reference),
  argmax agreement, DLR agreement, hostless kept/rejected at the fixed
  cuts, ECE.

Appendix A: universes (clean3, half-light DLR, ext_none, truth_nearest,
Legacy candidates, galaxy-only DP2) as one table of headline numbers.
Appendix B: the DC2 PA mirror and offset diagnostic. Appendix C: Tractor
and SExtractor conventions.

---

## 4. Inconsistencies found

Ordered by how much of the story they touch.

### 4.1 The headline universe is not the spec default (sample cut)
`baseline` runs `clean3`; the spec default and the only justified option
is `physical`. Every DP0.2/OU2024 headline figure and number quoted so far
(97.42 / 91.28, 74.4 / 71.7, the 7.30 % hostless, the ROC optima 0.0215)
is on clean3; the physical numbers (98.32 / 95.89, 75.1 / 74.6, hostless
2.35 %, J optimum 0.032) are the ones the rationale argues for. Fix: make
`physical` the baseline (Section 1). All physical outputs already exist.

### 4.2 Two DLR definitions, and the headline one is the self-referential one
`half_light` is the headline, `moments` the literature. On DP2 this makes
the central SBW-vs-DLR figure degenerate (AGENTS: "near-identity is by
construction"); on DES the faithful DLR is in a different universe
(`des_faithful`), so the text quotes DLR 97.10 % and 97.60 % for the same
sample. On DP0.2 half-light DLR is 2 points *better* than moments and on
OU2024 4 points worse, which will read as cherry-picking if the two tiers
are shown under different DLRs. Fix: moments everywhere in the paper; the
half-light DLR becomes an appendix line "a DLR on the fitted ellipse does
X". Legacy candidates have no moments, so the Legacy robustness paragraph
states it uses the fitted ellipse.

### 4.3 Denominators differ between figures of the same tier
Section 2. Concretely: the matrix is on cat1+cat2 survivors but the
hostless ROC's "real SNe kept" is on the full sample including SNe whose
host is absent from the catalogue, so part of the "SBW hostless" loss in
the ROC is cat1, not a method property. The miss accounting uses yet
another N. Fix: ROC real-SN curve on rung 4 too (the SNe for which a host
could have been found), or on rung 2 with the cat1+cat2 fraction drawn as
the ceiling; either way one choice, applied to all four tiers. The DES
matrix should be on survivors like the sim matrices.

### 4.4 Random-position designs differ by tier
DP0.2/OU2024: one position per truth SN, uniform over the tract, redrawn
until a 1 arcmin cell holds an object (386,248 for 306,332 SNe, so the
random N is not the sample N). DP2: ten per transient inside its own
patch. DES: ten per SN in a 45 to 90" annulus around the SN. These give
different "random rejected" baselines (DP0.2 85.7 %, DP2 83.0 %, DES
75.8 % at flb 0.01) that are then compared across tiers on Fig 5. Fix:
one design (ten per SN in an annulus around the SN is the most
field-matched and is what the real tiers use; the sims can afford it
cheaply), and the same N convention.

### 4.5 Two headline metrics for the same quantity
DP0.2 reports mean P_SBW(true) with argmax secondary; DES and DP2 report
argmax agreement first and mean P in passing; DP2-vs-DELIGHT reports
"same host" fractions only. Fix: every tier reports, in this order, mean
P_SBW(reference), the four P bands, argmax agreement; DELIGHT gets the
same P bands of the DELIGHT object (they exist in
`dp2_cross_method_metrics`).

### 4.6 The DP2 tier uses two samples, three candidate catalogues and five universes
rc_salt (690) for SBW/DLR, parsnip_ddf (7,823) for DELIGHT, then clean
(2,601) and galaxy (2,055), each with DP2 Object and Legacy candidates.
Hostless ROCs, matrices and galleries exist for all of them. Fix: the
paper shows rc_salt for Sections 5 and 6 and dp2_ddf_clean for Section 7,
with Legacy and galaxy as sentences. The remaining universes stay in the
spec for the appendix but produce no paper figure.

### 4.7 Thresholds: fixed in the spec, tuned in the text
The spec fixes flb 0.01 and d_DLR 4; the notes report per-tier optima
(flb 0.0215 / 0.010 / 0.091 / 0.133; d_DLR 2.7 / 3.8 / 2.2 / 2.1 to 2.7)
and the `flb02` option was added after DP0.2. The paper must choose: a
fixed operating point tested for transfer (recommended, and it is the
spec's rationale) and never a per-tier optimum in a headline number. The
`flb02` option can stay as a universe for the transfer argument.

### 4.8 DES search radius and the published-host matching
`des_faithful` uses r15 to mimic DES; nothing else in the paper does. With
moments DLR as the paper default, the only remaining difference is r15,
which should be dropped (the DLR 97.60 vs 97.10 difference is the
ellipse, not the radius; verify by materializing `des_*` under moments at
r30, i.e. the baseline-with-moments universe).

### 4.9 Extendedness on real data is a sample cut in disguise
On DP2 the ext05 cut on candidates is what turns the DDF ParSNIP sample
from "DELIGHT finds a host for 90 %" into "SBW finds one for 47 %". The
cleaned sample removes stars/AGN from the *transient* list, which is the
right place, but the paper must say that ext05 then still removes faint
compact galaxies under 634 of the remaining transients. One
sentence on `ext_none` for DP2 (not materialized) would close this; it is
the same question as the 18 % OU2024 loss and the DES star-classed hosts,
so it deserves one paragraph covering all three tiers.

### 4.10 Spec text and outputs that no longer match the project
- The analysis description promises five tiers including the DELIGHT
  visually matched set; it is unavailable. `delight_*` outputs and
  `delight_catalogue` have no recipe and never will here.
- `summary_table` has no recipe; it is the natural Table 2.
- `dp02_credible_set`, `dp02_pa_mirror`, `*_band_dependence`,
  `*_profile_sweep` are figures/tables that the paper will cite as
  numbers; fine to keep as outputs, but they should not be described as
  paper figures in the output descriptions.
- `candidate_quality` panel (c) of the miss-accounting figure is a
  universe sweep inside a figure; the convention elsewhere is one
  universe per figure.
- Vocabulary drift: `des_recovery`, `ou24_recovery` are agreement/recovery
  respectively; the DES figures rewrite "true" to "published" at save
  time. Keep the ids but make the captions use recovery (sims) /
  agreement (real) consistently, and consider renaming `des_recovery`.

### 4.11 The "DP0.2 optima" used for the transfer test are baseline's, in every universe
Every `*_hostless_metrics.json` carries `dp02_optima_transfer.dp02_optima`
= {flb 0.0215, d_DLR 2.7, credible t 0.034} and every credible-set JSON
`dp02_knee_transfer.t` = 0.034, including the physical, sep10 and
truth_nearest DP0.2 files whose own optima are 0.032 / 2.5 / 0.046,
0.018 / 2.9 / 0.028 and knee 0.042. The values are constants in the
recipes, not read from the DP0.2 output of the same universe. Once the
paper universe is physical, the transfer test on OU2024, DES and DP2 must
mark physical's optima (0.032, 2.5, 0.046), and the open diamonds on
every ROC move. Fix: pass the DP0.2 hostless/credible JSON of the same
universe as an input to the other tiers' metrics recipes.

### 4.12 One DES DLR rate disagrees with itself
In `des_recovery.json` the headline `dlr_recovery` is 0.9710 (1,576 of
1,623) while `per_band_profile["r/exponential"].dlr_agrees` is 0.9729;
in des_faithful 0.9760 vs 0.9784. The two are the same quantity under
two denominators or two matching rules (most likely the sweep counts
agreement among hosted pairs only). Find which and keep one.

### 4.13 Result files do not record the universe they were made under
The `settings` block of the recovery / miss / hostless JSONs lists band,
profile, radius, d_DLR, extendedness, quality and the two thresholds, but
not `dp02_sample_cuts`, `truth_match`, `dp2_sample`,
`dp2_candidate_catalogue` or `reference_match` (the last two appear only
in the cross-method file). So physical, sep10 and truth_nearest carry
settings blocks identical to baseline's. The manifests hold the
decisions, so provenance is intact, but a figure or JSON read on its own
cannot say which sample it is on. Fix: write every active decision into
`settings` from the recipe flags (they are all already CLI flags).

### 4.14 Real-data and simulated tiers report different quantities
Miss accounting, credible sets and the truth-match comparison exist only
for the simulations (correct, they need truth), but DES could carry the
same cat1/cat2 accounting against the published host (14 absent, 15
point source, 1,594 kept are already in `published_host_status`) and the
same P-band vocabulary, and today DES keys them `sbw_argmax_recovery` /
`mean_p_sbw_true` while DP2 keys them `same_host` / `mean_p_sbw_of_dlr_host`.
Confidence intervals come in three formats (`*_ci68` siblings, nested
`ci68`, `lo68`/`hi68`) and the DP2 disagreement and credible-set files
have none. A shared `summary_table` recipe is the forcing function: write
it first and make every tier emit the keys it needs.

### 4.15 Smaller things seen in the files
- Random-to-real ratios differ: 1.26 (DP0.2), 1.99 (OU2024), 10 (DES,
  DP2); see 4.4.
- OU2024's SBW Youden optimum sits exactly on the default grid point
  0.01, so "the DP0.2 cut transfers at a cost of 0.003" is a statement
  about a 121-point log grid; say the grid.
- `kang_crosscheck`: on DP2 Object candidates Kang's own nearest match
  finds 408 / 229 / 200 objects (ddf / clean / galaxy) that the 1"
  `reference_match` does not; on Legacy the two agree exactly. The 3"
  option matches 97.8 %. Worth one sentence, not a decision change.
- The median SN-to-DELIGHT separation rises from 0.18" (ddf) to 0.54"
  (clean) to 0.93" (galaxy). This is expected (the removed transients are
  the ones where DELIGHT marked the point source under the transient) and
  is a good one-line argument for the cleaning.
- `dp2_p_top_host` exists only as a figure; its numbers (median P_max
  1.000, 99.0 % above 0.9) should be in a JSON if the text quotes them.
- The miss-accounting definition says cat3 is `flb_best < 0.01` while the
  hostless rule says "no SBW host or flb < t"; make the wording the same
  (an SN with no surviving candidate is hostless under both).

### 4.16 Numbers quoted in notes that predate the current chain
AGENTS.md still carries the first-results paragraph (radius fallback,
97.42 vs later 97.42 by coincidence), the prototype DES 98.6 / 98.4, and
OU2024 headline numbers computed before the no-fallback change. None
should reach the paper; only numbers read from current result files do.

---

## 5. Decision audit summary

Justified and kept as is: `sbw_profile`, `band`, `search_radius` (r30),
`dlr_max`, `extendedness_cut`, `sbw_hostless_cut` (as a fixed point),
`p_true_threshold`, `truth_match`, `reference_match`, `dp2_sample` (two
samples, one reason each), `dp2_candidate_catalogue`.

Justified but the wrong option is in the headline: `dp02_sample_cuts`
(use physical), `dlr_radius` (use moments).

Not a paper decision: `candidate_quality` (keep none; one sentence),
`exclude_suspect` (never varied), `delight_catalogue` (tier removed).

Options never materialized and not needed: `dlr_max: d5`,
`search_radius: r60`, `band: g/i/z` as universes (the band sweep is an
output), `reference_match: pos_3arcsec / frac_semimajor` (one sentence
from the DELIGHT metrics), `dp2_sample: md_final / rc_all`.

---

## 6. Work list to reach the consistent set

1. Spec: done 30 Oct 2026 (moments default, baseline = spec defaults,
   appendix universes under their old ids, DELIGHT tier removed, validated,
   uncommitted).
2. Random positions: one design for all tiers (ten per SN in a field
   annulus); rerun the random chains (DP0.2/OU2024 ~25 min per universe on
   NERSC; DP2 minutes; DES minutes on the laptop).
3. Figures: hostless ROC real-SN curve on the chosen rung; DES matrix on
   survivors (`--survivors-only` already exists in `dp02_fig.py`, DES runs
   through it); drop panel (c) of the miss figure; a four-tier ROC figure
   and a three-way DELIGHT matrix figure as new outputs; `summary_table`
   recipe.
4. Materialize baseline (new meaning) for all four tiers; DES on the
   laptop, the rest on NERSC under the usual allocation rules.
5. Metrics plumbing: the DP0.2 optima read from the same universe's
   output (4.11); every decision written into `settings` (4.13); one CI
   format and one key vocabulary across tiers, driven by the
   `summary_table` recipe (4.14); resolve the DES DLR pair (4.12).
6. Rewrite the headline paragraphs of AGENTS.md from the new outputs and
   delete the superseded numbers.
7. Then `index.md`.
