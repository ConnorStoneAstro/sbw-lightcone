# SBW host association performance — Lightcone project

This directory is the Lightcone/ASTRA analysis behind a paper on the
surface-brightness-weighted (SBW) host-association method for transients.
`astra.yaml` is the source of truth for inputs, outputs, decisions and the
literature evidence; `universes/baseline.yaml` selects the package defaults;
`index.md` is the MyST report. The method itself lives in the sibling package
repo `../hostassoc-repo` (pip name `sbw-host`, package `sbw_host`, GitHub
`ConnorStoneAstro/sbw-host`); this project only consumes it.

## Day to day

```bash
uvx astra-tools@0.2.18 validate                     # after every spec edit
uvx astra-tools@0.2.18 validate astra.yaml --verify-evidence   # quotes vs cached PDFs
lc status --json                                    # what is current / stale / no recipe
lc compute launch --wait                            # local 8-CPU cluster named "local"
lc run local -- python src/<script>.py ...          # probe a script in the sandbox
lc materialize local <output_id>                    # make one output for real
uv add <pkg>                                        # the only way to change the environment
```

Commit by path (`git add src/ astra.yaml AGENTS.md`), never `git add -A`.

## Project Notes

Scoping was done on 2026-10-26 with Connor Stone. State of play:

- **Framing.** SBW is the assumed framework, not a contender. The paper asks
  (a) which *form* of SBW works best (exponential fit, free-n Sersic fit, or
  deblended image stamps, all in `sbw_host`), and (b) how SBW relates to DLR:
  disagreement rate, band dependence, agreement versus d_DLR. "Which method
  wins" is not the question. A cheap cosmology propagation (wrong host, wrong
  redshift, delta mu) is wanted if it stays cheap; the full SNANA/BBC route is
  out of scope. Qu+24 (delta w = 0.0013 +/- 0.0026 for < 2% mismatch) is the
  benchmark for that.
- **Form-of-SBW is both a decision and an output.** `sbw_profile` sets the
  headline form for every tier; `dp02_profile_sweep` computes all three forms
  side by side on the truth sample and does not depend on that decision.
  Band dependence likewise is an output (`dp02_band_dependence`,
  `dp2_band_consistency`), with `band` only fixing the headline band.
- **Simulation tiers: DP0.2, and OpenUniverse2024 is coming back (27 Oct 2026).**
  OU2024 was dropped on 26 Oct on the belief that its WFD coadds had been
  purged; that was wrong. The roman-desc-sims butler repo on NERSC
  (`/global/cfs/cdirs/lsst/production/gen3/roman-desc-sims/repo`) holds,
  besides the multiband Object table (`u/descdm/step3b_wfd_000..066_w_2024_22`,
  26 tracts, 70 GB; CModel bulge+disk / moments / Kron, no MultiProFit
  Sersic), `deepCoadd_calexp` for 22 tracts (683 GB) in
  `u/descdm/step3_wfd_000..062` and `step3a_wfd_063..066`, plus scarlet
  model data and DIA tables, so the stamp form of SBW is possible there. The
  user decided OU2024 is folded back in as a second simulated tier; it is
  wired on branch `ou24-tier` (next item).
- **OU2024 tier (27 Oct 2026, branch `ou24-tier`, NERSC clone
  `$PROJECT/host_association/paper-ou24`; merged into main as 26031fe, run from the main clone since 28 Oct).** Inputs are
  three index csvs under `data/`: `ou24_object` (67 objectTable files, all
  needed because each tract's runs hold disjoint patches),
  `ou24_sn_truth` (33 skyCatalog `snana_*.parquet`, the SNANA dump
  `OU24_TRANSIENT_TRUTH_HEAD.DAT.gz`, 155 `*_HEAD.FITS.gz`) and
  `ou24_host_truth` (skyCatalog `galaxy_*`, `galaxy_flux_*`,
  `galaxy_sed_*.hdf5` for 33 pixels + the throughputs), written by
  `src/ou24_index_objects.py` / `src/ou24_index_files.py`. Chain:
  `ou24_sn_table` (521,766 transients, all hosted, 5 min; B/T from the
  component SEDs) -> `ou24_truth_match` (hand-built; `truth_match:
  matches_truth` = magnitude-aware, `nearest_1arcsec` = nearest; 51.7% / 61.1%
  of hosts) -> `ou24_sample` (same `dp02_sample_cuts`; FIXMAG always dropped;
  all z kept, median 1.5) and `ou24_candidates` (from the whole sn_table) ->
  the DP0.2 scripts unchanged (`--tier ou24` only changes titles;
  `dp02_random_positions.py --sne` draws over a tract's file union;
  `dp02_offsets.py` without `--shapes` reads the sample's offset columns).
  No `--radius-fallback` anywhere on OU2024. Whole chain ~45 min per universe
  on 16 CPUs. New: `ou24_recovery_by_host_mag`, `ou24_dp02_comparison`
  (takes `dp02_associations`, which already carries host_mag_r, so
  `dp02_sample` is not an input), `ou24_pa_check(_metrics)`. Materialized:
  baseline, physical, and the new universe `ext_none` (OU2024 outputs only;
  its DP0.2 outputs are unmade, and `ou24_dp02_comparison` was left out
  there because it would pull the DP0.2 chain).
  * *Headline (baseline, r, exponential, 262,710 transients):* SBW 74.4%,
    DLR 71.7%, agreement 70.4%, DLR hostless 24.3%, true host among
    candidates only 81.7%, ECE 0.24. z < 1: SBW 85.7%, DLR 80.6%. Hosts
    r < 24 (58k): SBW 98.5%, DLR 93.2%, agreement 93.9% (z < 1: 98.7 / 93.0),
    against DP0.2 97.9 / 91.3 at r < 24. Reweighted to DP0.2's host-r
    distribution OU2024 gives SBW 97.2%, DLR 92.4% (DP0.2 97.4 / 91.2). So
    the machinery agrees; the headline is low because OU2024 hosts are
    faint (median sample host r 25.1 vs 22.3, the host-mass bug) and because
    at r 24-26 OU2024 recovers much worse than DP0.2 at the same magnitude
    (SBW 0.89 at 24-25, 0.65 at 25-26, 0.40 at 26-27; DLR > SBW past r 26).
    Likely causes: ext05 (below) and OU2024's 15 s exposures (coadds ~0.4
    mag shallower than intended; insight `ou2024_lsst_exposures_15s_shallower`).
  * *physical (235,637):* SBW 75.1%, DLR 74.6%; r < 24 SBW 99.0%, DLR 96.7%
    (DP0.2 physical 98.4 / 95.9); reweighted SBW 97.4%, DLR 95.1%.
  * *Extendedness.* ext05 drops the true host for 18.3% of transients
    (faint compact hosts with refExtendedness < 0.5). `ext_none`: true host
    always a candidate, SBW 85.2% (+10.8), DLR 81.2% (+9.5), z < 1 92.3 /
    86.3, r < 24 unchanged (98.8 / 92.9), ECE 0.13; the cost is random
    rejection SBW 91.8 -> 90.3% and DLR 86.8 -> 83.8% (flb 0.01 / d 4). The
    gain is all in faint hosts; kept as a universe, not the default.
  * *Hostless (baseline):* random positions rejected SBW 91.8%, DLR 86.8%
    (2x2: both reject 86.6%, SBW alone 5.1%, DLR alone 0.2%); real kept at
    flb 0.01 70.8% vs DLR d4 75.7%, d3 70.8% / 91.2% rejected, so SBW flb
    0.01 and DLR d3 sit at nearly the same point, as on DP0.2. 99.3% of
    random positions still get a P > 0.5 host. Real SBW hostless 29% is
    mostly the faint hosts whose true host is not a candidate.
  * *Matrix (baseline, flb 0.01):* SBW right / DLR wrong 2,207, DLR right /
    SBW wrong 399, both wrong 4,427 (4,139 same object), SBW hostless / DLR
    true 10,784.
  * *Profiles (r):* baseline exponential 74.4%, de Vaucouleurs 74.7% (ECE
    0.14), moment Gaussian 76.3%; ext_none 85.2 / 87.4 (ECE 0.04) / 91.5%.
    Unlike DP0.2 (de Vaucouleurs best), the moment Gaussian wins on OU2024,
    because it keeps faint hosts whose CModel fit is poor.
  * *Bands (baseline):* SBW u 31%, g 64%, r 74%, i 66%, z 55%, y 47%, DLR
    ~0.65-0.70 flat; the SBW host differs from r for 63% in u. Cause: the
    true host's light model is invalid (no usable HSM moments in that band)
    for 61% in u, 18% g, 7% r, 18% i, 33% z, 43% y, so it drops out of SBW;
    for hosts r < 24 it is valid in g-y >= 98.6% (u 76%). A shape from the
    reference band when a band's moments fail would fix it; not done (a
    code decision for Connor, cf. the no-fallback rule for radii).
  * *PA check: no mirror bug in OU2024.* For q < 0.5, offset > 1 R_e
    (93,864): <cos 2(psi - beta)> = 0.306 +/- 0.002, mirror 0.000; against
    the pipeline moment axis of the matched object (45,660) <cos 2(psi -
    phi)> = 0.177, mirror -0.004; beta vs phi 0.48. beta = 0.5 atan2(g2, g1)
    = SNANA a_rot is measured from +x West toward North in the sky_to_plane
    frame, i.e. the same convention as the pipeline moments (the prep note
    "beta = -phi" came from an East-based frame). `ou24_sample` offsets use
    +beta.
  * *Offsets:* failure rates flat (~22-25%, set by faint hosts) to ~2 R_e of
    the dominant component, then rise; both-wrong-same ~2% rising to 6-10%
    beyond 5 R_e. The physical cut removes 10% of the clean3 sample.
  * *Operational.* `ou24_dp02_comparison` uses the DP0.2 results as of
    acb147c (with radius fallback); remake after merging main. The NERSC
    clone has remote `mainclone` = `../paper`, used once for `git annex get`
    of the DP0.2 associations. An ssh disconnect detaches the `lc` client
    but the login-node process keeps going.
- **DP0.2 runs entirely on NERSC from parquet.** Inputs: DESC merged truth
  `DC2-prod/Run2.2i/truth/truth_merged_summary_v1-0-0/match_dr6_v2/truth_tract*.parquet`
  (166 tracts; SNe = truth_type 3, 445k; host rows = truth_type 1 with mag_r;
  its match_objectId is the DESC DR6 processing and must not be used),
  Rubin's DP0.2 truth match `shared/rubin/DP0.2/match_ref_truth_summary_objectTable`
  (id string to DP0.2 objectId; 93% of SN hosts in tract 2897 have one), and
  `shared/rubin/DP0.2/objectTable` (157 tracts, 1 TB). TruthSummary itself is
  not on NERSC (TAP only), which is why the DESC merged truth is used. No
  DP0.2 images exist on NERSC, so the stamp form stays parked (RSP cutouts
  on a subsample later).
- **DP0.2 column conventions (src/dp02_lib.py).** `bdReB`/`bdReD` are
  determinant radii `(xx*yy-xy^2)^(1/4)` = sqrt(a b) of the de Vaucouleurs /
  exponential half-light ellipse in 0.2" pixels; semi-major = r/sqrt(q). The
  v23 `E1` functor had no parentheses (`xx - yy/(xx+yy)`), so `bdE1` is
  garbage and q, phi come from PSF-subtracted HSM moments
  (`ixx-ixxPSF` ...), falling back to raw moments. No fracDev is exported,
  so there is no bulge+disk composite: on DP0.2 `sbw_profile` maps
  exponential -> n=1 with bdReD and the total cModelFlux, sersic -> n=4 with
  bdReB and cModelFlux, moments -> Gaussian (n=0.5) with the moment
  ellipse. DLR `half_light` uses the bdFlux-weighted mean of bdReB/bdReD;
  `moments` uses the rms semi-axes. Checked on synthetic galaxies (m, d_DLR
  on major/minor axis) before running.
- **DP0.2 pipeline shape.** `dp02_candidates` (one row per SN-object pair
  within the search radius, 16 workers over the tract files) feeds
  `dp02_sample` (truth + host objectId + cuts) and `dp02_associations`; the
  per-band and per-profile sweeps (`src/dp02_sweep.py`) recompute the
  association in-process and do not depend on the headline band/profile.
  `dp02_hostless_indicator` still has no recipe (needs random-position
  candidates). Figures asked for on 26 Oct: SBW reliability curve
  (fraction of top hosts that are true vs P_SBW, 1:1 if calibrated) per
  band and per profile, the SBW/DLR/truth agreement matrix, agreement vs
  d_DLR; added on my side: recovery vs redshift and the P(true host)
  distribution. `lc status` hashes the declared input directories by
  content, so the first status after declaring the 1 TB objectTable is slow.
- **DP0.2 quirks** (from the 22 Oct 2026 prototype): ~5% of SN truth rows point
  to an empty host row; ~9% sit > 3 arcsec from the declared host; all z > 1
  and 10% of z < 1 SNe are deliberately hostless (control sample, not
  failures). Keep objectIds int64 end to end; sync TAP joins time out, use
  chunked `id IN (...)` queries or async.
- **First DP0.2 results (26 Oct 2026, baseline universe, 306,332 clean SNe,
  all z < 1 because DC2's z > 1 SNe are hostless).** SBW top host = true host
  97.42% (r, exponential), DLR 91.24%, agreement 92.62%, DLR hostless 3.9%,
  true host among candidates 99.7%. Profiles: de Vaucouleurs 97.76% (best,
  ECE 0.013), exponential 97.42% (ECE 0.015), moment Gaussian 96.87% (ECE
  0.025). Bands: u 92.6% (ECE 0.06), g 96.9%, r/i/z/y 97.4-97.7%; the SBW
  host changes vs r for 7% of SNe in u and < 1% elsewhere. SBW is
  *overconfident* in the 0.5-0.95 range (reliability curve ~0.1 below 1:1)
  and converges at P -> 1; ECE is driven by that band. Agreement matrix:
  SBW right / DLR wrong 3.1%, SBW right / DLR none 3.3%, DLR right / SBW
  wrong 0.2%, both wrong 1.7% of which 4,851/5,305 pick the *same* wrong
  object (likely the truth match pointing at a different deblended child;
  the `truth_match: nearest_1arcsec` universe would test this). Recovery
  falls with separation (SBW 0.73 at 3", DLR 0.37) and at both magnitude
  ends (bright 16-18 mag: DLR ~0.8, SBW ~0.95; faint > 25). Agreement
  with DLR declines with the SBW host's d_DLR and collapses past 4 by
  construction. Run cost on `shared-interactive`: candidates 9 min on 16
  CPUs, the whole 16-output chain 42 min.
- **DP0.2 hostless indicator (27 Oct 2026, baseline).** `dp02_random_positions`
  draws one position per truth SN uniformly (RA, sin dec) inside its tract,
  redrawn unless a 1 arcmin cell holds a primary DP0.2 object (386,248
  positions in 152 tracts; 7,774 SNe sit in the 13 truth tracts with no
  objectTable file and get none); ids are the negated truth ids. Random
  candidates (40.5M pairs, 10 min) and associations go through the same
  scripts (`--positions`). `flb_max` (max frac_light_beyond over candidates)
  was added to `associate_one`; `dp02_associations` was rematerialized for it
  (recipe gained `--workers 16`), headline numbers unchanged. Results: SBW
  P_SBW is a *relative* probability, so 98.9% of random positions still get a
  top host with P > 0.5 (79% with P > 0.9); hostlessness must come from
  flb. At flb_best < 0.01: 92.7% of real SNe kept, 85.7% of random rejected,
  14.2% of random kept on a P > 0.5 host; DLR d <= 4: 96.1% kept, 76.3%
  rejected (23.5% random kept with P > 0.5); d <= 3: 92.0% / 84.8%. On the
  ROC SBW lies slightly above DLR everywhere (about 1 point of random
  rejection at equal real kept); flb_max is indistinguishable from flb_best.
  Expected: for one exponential profile flb at 4 r_e is 0.0094, so flb 0.01
  and d_DLR 4 are nearly the same ellipse; the operating points differ
  mostly because DLR's radius is the bulge+disk flux-weighted one. Half the
  random positions have flb_best < 1e-5. 3x3 matrix (flb 0.01): SBW hostless
  7.30% of real SNe, of which 2.42% DLR-true, 1.05% DLR-wrong, 3.82% DLR
  hostless; SBW right 88.71% overall (2.35% with DLR wrong, 0.08% DLR none);
  recovery with the cut 91.14%. `dp02_hostless_indicator.pdf`,
  `dp02_hostless_metrics.json`. On NERSC the clone's hub remote is `origin`
  (`nersc` only on the laptop); `lc compute launch` there gives a 1 h
  allocation by default, enough for this chain (~25 min).
- **Truth match, sample cuts, Sersic tails, cases (27 Oct 2026).** All numbers
  below are after the `--radius-fallback` change (next item); r band,
  exponential, flb 0.01, d_DLR 4.
  * *Truth matchers* (`truth_nearest` universe, `dp02_truth_match_comparison`).
    `dp02_sample` now always carries `host_objectId_mt` and `host_objectId_nn`
    (+ `match_sep_arcsec`, median 0.009"); the parquet metadata key
    `truth_match_precut` counts them before the matcher cut. Of 321,161 hosted
    SNe within 3" of a host row: both match 306,332, MatchesTruth-only 0,
    nearest-only 4,464, neither 10,365, and they name different objects for
    only 75 (70 of those MT objects are other primary candidates, 5 are
    non-primary). Where they agree (306,257) SBW 97.54%, DLR 91.25%. The
    same-wrong-host cell is *not* a truth-labelling artefact: on the common
    set it is 4,177 against MatchesTruth and 4,165 against nearest.
    `truth_nearest` overall (310,796 SNe): SBW 97.40%, DLR 91.11%,
    both-wrong-same 4,275.
  * *Sample cut options* are `clean3` (old `clean`), `clean10`, `physical`
    (spec default) and `all`; `baseline`, `truth_nearest` and `dlr_moments` stay
    on `clean3`. `sep10` (333,265 SNe): SBW 93.02%, DLR 85.07%, agreement
    88.55%, ECE 0.053, both-wrong-same 10,437. The 26,933 SNe at 3-10" are
    recovered ~42% by SBW and ~15% by DLR.
  * *cosmoDC2 shapes.* `cosmodc2_hosts` is an index csv
    (`data/cosmodc2_v1.1.4_image_index.csv`, written by
    `src/dp02_index_objects.py <dir> <out> <glob>`) of the 131 healpix files in
    `/global/cfs/cdirs/lsst/shared/xgal/cosmoDC2/cosmoDC2_v1.1.4_parquet_r2`
    (2.9 TB; `cosmoDC2_v1.1.4_parquet` is a symlink to it). truth
    `host_galaxy` = cosmoDC2 `galaxy_id`; 354,095 of 394,021 SN host ids found
    (the rest are the empty host rows); the extraction takes 1-4 min.
    79% of the SN hosts have B/T_i > 0.5 (stellar-mass-weighted host choice).
  * *DC2 position angles are mirrored.* DP0.2 HSM moments of large elongated
    hosts give PA (East of North) = -`position_angle_true` (median |dPA|
    0.8 deg over 533 hosts), but SN offsets align with +PA
    (<cos 2(psi-PA)> = 0.59 for q < 0.3, 0.005 for -PA). SN radii follow the
    semi-major R_e exactly (median r/R_e 1.00) on a rounder ellipse (axis
    ratio ~ q^0.4; KS 0.006 for pure bulges). So DC2 placed SNe on the mirror
    image of the rendered galaxy. It costs little: flat hosts fail *less*
    (both-wrong-same 1.4% for q < 0.3 vs 1.6% for round). `host_offsets` in
    `dp02_lib.py` uses +PA (the placement convention). Recorded in the
    docstring and the dp02_sample_cuts rationale, not yet as a prior insight
    (no literature quote).
  * *Offset diagnostic* (`dp02_offset_diagnostic`, sep10 / baseline / physical):
    every failure rate is flat to ~1.5 R_e of the dominant component, then
    climbs. Bulge-dominated, sep10: both-wrong-same 0.5% at 1.1 R_e, 2.7% at
    2.0, 6.3% at 2.7, 12% at 3.6, 18% at 4.8, ~25% plateau beyond 6; DLR
    wrong/none 16% -> 38% -> 75% over 2.7-4.8 R_e. Disk-dominated rises
    later and lower (4% at 3.6). In the baseline, 84% of same-wrong-host SNe
    lie beyond 2 R_e, 65% beyond 3. So the cell is DC2's untruncated n = 4
    bulge tails, and the cutout gallery shows those SNe sitting on a
    neighbouring galaxy. The DC2 paper (IOP DOI 10.3847/1538-4365/abd62c, now
    cached) normalises offsets by the hosting component's semi-axes and
    mentions no truncation: insight
    `dc2_sn_offsets_scale_with_hosting_component_size`.
  * *Physical cut* (`N_RE = 3.0` in `src/dp02_sample.py`, bulge R_e for
    B/T_i > 0.5 else disk R_e, elliptical with +PA; hosts without shapes are
    kept, none in practice): 66,269 SNe beyond 3 R_e (18.7% of the 354k with a
    host row). `physical` (274,633 SNe): SBW 98.41%, DLR 95.86%, agreement
    96.52%, ECE 0.008, SBW hostless 2.22%, DLR hostless 0.57%,
    both-wrong-same 2,234; recovery stays ~0.6-0.7 for SNe 5-10" from large
    hosts. N = 4 or 5 would remove ~12% / ~8.5% of the sep10 sample.
  * *Radius fallback: removed (27 Oct 2026, later).* It briefly fell back to
    the other CModel radius, then the moment radius, for hosts with NaN
    `r_bdReD`. Now a candidate with no usable radius or flux for the profile
    is *always dropped* before association, and a true host dropped that way
    counts as cat2 "unusable" in the miss accounting. All numbers in the
    first-results paragraph above predate this; current ones are in *Miss
    accounting, quality cuts, ...* below.
  * *Agreement matrix with random positions* (peer request):
    `dp02_agreement_matrix` and `dp02_recovery` take `dp02_random_associations`;
    the figure adds a 2x2 for random positions and `dp02_recovery.json` gains
    `agreement_matrix` and `random_matrix`. Random positions at flb 0.01 /
    d_DLR 4 half-light: 75.8% rejected by both, SBW alone 9.8%, DLR alone
    0.5%, neither 13.9% (SBW rejects 85.7%, DLR 76.3%). Every universe that
    has a matrix now also computes its own random chain (~25 min each, since
    per-universe outputs are never shared).
  * *DLR on the moment ellipse* (`dlr_moments` universe, peer request): DLR
    89.28% (half-light 91.24%), DLR hostless 8.00% (3.91%), agreement
    90.44%. The matrix moves DLR failures from "wrong" to "hostless": SBW
    right / DLR hostless 3,710 (half-light 444), SBW right / DLR wrong 3,322
    (7,340). Hostless ROC: moments d <= 3/4/5 keep 85.5/92.0/95.9% of real
    SNe and reject 91.2/84.9/77.4% of random positions, against half-light
    92.0/96.1/98.0% and 84.8/76.3/67.2%, and SBW flb 0.01 at 92.8% / 85.7%.
    So moment-DLR at d 4 sits at nearly the same operating point as SBW flb
    0.01, slightly below it, and all three lie on almost the same ROC curve.
    The random matrix there: SBW alone rejects 3.9%, DLR alone 3.2%.
  * *Cases and gallery.* `dp02_disagreement_cases` (baseline; seed 20261027,
    12 per class) needs only one universe because the sample carries both
    matchers. Classes (baseline, p_true_threshold 0.05): truth_matchers_disagree
    75, sbw_excludes_true_host 2,467 (P(true) < t), dlr_wrong_or_hostless
    26,726, sbw_most_probable_ne_dlr 22,239, true_host_cut_by_object_table
    1,873 (cat2). `src/dp02_gallery.py` began as a copy of
    `../dp02_cutouts_proto/gallery.py` (sha256 e6c7f5a3...880b2); only its
    labels changed since. The gallery now *materializes on NERSC*:
    `rsp_token` is `/global/cfs/cdirs/m1727/cstone/.rsp_image_tap_token`
    (mode 600; never print it), compute nodes reach data.lsst.cloud, and the
    recipe caches cutouts under `$TMPDIR`.
  * *Operational.* The NERSC clone must `git fetch origin && git annex merge`
    before `git push origin main git-annex`, or the git-annex push is
    rejected. A remade output whose bytes are identical leaves its downstream
    only "behind", not stale. `lc compute launch ... --time 4h` gives the full
    shared_interactive window. The full five-universe rematerialization with
    random chains took 1 h 50 min on 16 CPUs. `$PROJECT` paths inside
    `n.sh`-style double-quoted ssh commands expand on the laptop, so quote
    remote globs in single quotes.
- **Miss accounting, quality cuts, P(true), credible sets, Youden J, PA
  mirror (27 Oct 2026, evening; no radius fallback, candidate_quality none,
  p_true_threshold p005, r, exponential, flb 0.01, d_DLR 4 half-light).**
  * *Language.* SBW is scored by P_SBW(true host), not argmax. No figure says
    "SBW host = true"; the wording is "SBW most probable host". The argmax
    recovery is kept as a secondary number (`sbw_argmax_recovery`).
  * *Headline (baseline, 306,332 SNe).* Mean P_SBW(true) 0.9705 +/- 0.0003
    (argmax 97.42%, DLR 91.28%, agreement 92.74%). Bands: P > 0.5 91.13%,
    0.05-0.5 0.87%, < 0.05 0.70%, SBW hostless (flb < 0.01) 7.30%. Physical
    (274,633): mean P 0.9799, argmax 98.32%, DLR 95.89%, SBW hostless 2.35%.
    The 4x3 matrix: P > 0.5 / DLR wrong 6,919, P < 0.05 / DLR wrong 2,108,
    SBW hostless / DLR true 7,271.
  * *Miss accounting* (`dp02_miss_accounting`). The denominator is 321,161
    hosted SNe: the 14,829 whose true host has no DP0.2 object count as cat1.
    Baseline SBW: cat1 absent 4.62%, cat2 cut 0.58% (point source 989,
    unusable radius/flux 884), cat3 genuine 7.04% (P < t 1,977, hostless
    20,641), success 87.76%. DLR: cat2 0.58% (+884 true hosts with no DLR
    ellipse), cat3 7.74% (other closer 14,021, beyond d 4 10,827), success
    87.06%. Physical: SBW cat3 1.99%, success 92.85%; DLR 3.33% / 91.51%.
    Most of SBW's cat3 is the hostless cut, not a wrong pick.
  * *candidate_quality.* `basic` (cModel_flag or S/N < 5) removes 7.9% of
    candidates but only 277 true hosts (cat2 0.58 -> 0.67%); mean P 0.9705 ->
    0.9702, so it changes nothing. `strict` (+ blendedness > 0.5 or > 50
    siblings) removes 13.6% of candidates and 32,012 true hosts (cat2 10.6%),
    mean P 0.874. It lowers SBW cat3 only 7.04 -> 6.44% and DLR 7.74 ->
    5.39%. So object-table cuts are not a route to fewer failures on DP0.2:
    the true hosts are as blended as the impostors. Sibling count barely
    discriminates: 96% of true hosts are deblended children, median 10
    siblings vs 12 for wrong picks (`SIBLING_MAX = 50` removes 0.85% / 2.4%).
  * *Credible sets* (`dp02_credible_set_metrics`, set = candidates with
    P > t). Baseline at t 0.01: coverage 98.92%, mean size 1.087. At t 0.05:
    98.58% / 1.046. Knee at t 0.034 (98.68%, size 1.055). Top-k: k = 1
    97.42%, k = 2 99.15%, k = 3 99.32%. The ceiling is the 99.39% of true
    hosts that survive the cuts. Physical knee t 0.046 (99.04%, 1.035). So a
    one-extra-candidate set buys ~1.2 points over argmax at a mean cost of
    0.05 candidates.
  * *Youden J* (hostless indicator, real kept + random rejected - 1). On
    baseline the flb_best maximum is J 0.789 at flb 0.0215 (89.9% real kept,
    89.1% random rejected) and DLR's is J 0.778 at d_DLR 2.7 (89.9% / 87.9%).
    Physical: flb 0.032 (J 0.853), d 2.5 (J 0.849). flb_max matches
    flb_best. 0.0215 is > 2x from 0.01 and from 0.05, so `sbw_hostless_cut`
    gained `flb02` (0.02); the default stays flb01. Do not retune these on
    DP0.2 before OU2024; they are meant to transfer unchanged.
  * *Other universes.* truth_nearest (310,796): mean P 0.9693, argmax 97.30%,
    DLR 91.14%, SBW cat3 7.14%, knee t 0.042, J flb 0.0215 / d 2.7.
    sep10 (333,265): mean P 0.9260, argmax 92.94%, DLR 85.12%, SBW cat3
    13.56% vs DLR 13.64%, knee t 0.028 (cov 95.34%, size 1.093), J flb 0.018
    / d 2.9. dlr_moments: SBW as baseline, DLR 89.19%, DLR cat3 9.73%,
    J-optimal d_DLR 3.6 (J 0.773). The J-optimal flb is 0.018-0.032 in every
    universe, which is the case for flb02. Whole run: baseline+physical 2 h
    25 min, the other three 3 h 50 min on 16 CPUs (two 4 h allocations).
  * *PA mirror, as a prior-insight-like note* (`dp02_pa_mirror`, no
    literature quote yet). 7,664 elongated hosts (q < 0.5 in truth and in
    moments, SN beyond 1 R_e). SN offset angle against the pipeline moment
    PA: <cos 2 delta> = -0.002 (flat; median |delta| 45.2 deg). Against the
    mirrored moment PA (-PA): 0.326 (27.1 deg). Against the cosmoDC2 PA:
    0.335 (26.4 deg). median |PA_mom + PA_cosmoDC2| = 1.3 deg, while
    median |PA_mom - PA_cosmoDC2| = 44.4 deg. The cutouts show six large
    hosts with the SN on the mirrored ellipse, off the visible disk. So DC2
    rendered galaxies at -PA and placed SNe at +PA. Any offset-shape
    analysis on DC2 must use +PA (`host_offsets` does). The SBW/DLR
    association is unaffected because it uses the rendered image.
- **Aligned suite: four sequential miss categories, OU2024 on the DP0.2
  framing, PA radius check (28 Oct 2026; r, exponential, flb 0.01, d_DLR 4
  half-light, t 0.05, candidate_quality none).**
  * *Categories* (`dp02_miss.py --tier dp02|ou24`, `{dp02,ou24}_miss_accounting(_fig)`).
    Each category is conditional on passing the previous one: cat1 absent
    from the catalogue, cat2 object-table cut (point source / unusable
    radius or flux), cat3 no-host cut (SBW flb < cut; DLR nothing within
    d_DLR 4), cat4 genuine failure (SBW P(true) < t; DLR another candidate
    closer). Numbers are % of all hosted SNe (conditional % in brackets).
    DP0.2 baseline (321,161): absent 4.62, cut 0.58 (0.61); SBW hostless
    6.43 (6.78), genuine 0.62 (0.70), success 87.76; DLR hostless 3.37
    (3.56), genuine 4.37 (4.78), success 87.06. DP0.2 physical (287,784):
    SBW 1.70 / 0.29 / 92.85, DLR 0.17 / 3.16 / 91.51. OU2024 baseline
    (507,988 hosted, matcher-unmatched hosts count as absent): absent 48.28,
    cut 11.74 (22.69); SBW hostless 4.00 (10.01), genuine 0.16 (0.45),
    success 35.82; DLR 3.80 (9.50) / 0.86 (2.39) / 35.32. OU2024 physical
    (453,584): absent 48.05, cut 11.74; SBW 2.69 / 0.09 / 37.43, DLR 2.60 /
    0.67 / 36.94. So on both tiers SBW's losses are almost all the no-host
    cut and DLR's are wrong picks; OU2024's low headline is cat1 + cat2
    (faint hosts unmatched or cut), not association. `strict` quality on
    OU2024 lifts cat2 only 11.74 -> 13.87 (DP0.2 0.58 -> 10.55).
    OU2024 ext_none: cut 3.57 (6.90); SBW hostless 7.39, genuine 0.40,
    success 40.36; DLR 7.01 / 1.70 / 39.43; mean P 0.8491, knee t 0.042,
    J optima flb 0.012 / d 3.6 (DP0.2 flb 0.0215 -0.002, d 2.7 -0.015).
  * *P(true) and bands* (OU2024 now has them: `ou24_recovery`,
    `ou24_agreement_matrix --p-threshold`). Mean P_SBW(true): DP0.2 0.9705
    (physical 0.9799), OU2024 0.7419 (0.7494). Bands P > 0.5 / t-0.5 / < t /
    SBW hostless: DP0.2 91.13 / 0.87 / 0.70 / 7.30, OU2024 68.96 / 0.29 /
    1.59 / 29.16 (physical 71.83 / 0.22 / 1.43 / 26.53). `ou24_dp02_comparison`
    now plots SBW mean P(true) (remade against the current, no-fallback
    dp02_associations).
  * *Credible sets.* Knee t: DP0.2 0.034 (cov 98.68%, size 1.055), physical
    0.046; OU2024 0.051 (cov 74.97%, 1.080), physical 0.051. The DP0.2 knee
    0.034 on OU2024 gives coverage 75.02% / size 1.094: it transfers (the
    curve is flat between 0.03 and 0.05 on both tiers).
  * *Youden J and transfer* (`hostless_metrics` key `dp02_optima_transfer`;
    the ROC panel marks the DP0.2 optima as open diamonds). OU2024 baseline
    optimum flb 0.010 (J 0.626), d_DLR 3.8 (J 0.600); at DP0.2's flb 0.0215
    J is 0.003 lower, at d 2.7 0.017 lower, at d 4 0.001 lower. Physical:
    flb 0.018 (J 0.654; flb 0.01 and 0.0215 within 0.001), d 3.4 (d 2.7
    -0.008, d 4 -0.004). dlr_moments OU2024: d 3.8 (J 0.661), d 2.7 -0.016.
    Verdict: the SBW flb cut transfers anywhere in 0.01-0.02; DP0.2's DLR
    optimum d 2.7 does not, and the conventional d 4 is the safer default.
    Random rejection OU2024: SBW 91.8%, DLR half-light 87.8%, moments 91.6%.
  * *DLR half_light vs moments.* DP0.2: DLR 91.28 vs 89.19%, DLR hostless
    cat3 3.37 vs 7.25, genuine 4.37 vs 2.47. OU2024: 68.30 vs 72.57%, cat3
    3.80 vs 2.10, genuine 0.86 vs 0.35, success 35.32 vs 37.53. Moments
    help on OU2024 (no bulge+disk radius for faint hosts) and hurt on DP0.2.
  * *PA mirror radii* (`dp02_pa_mirror_metrics`, page 2 of `dp02_pa_mirror`).
    29,204 baseline hosts with q < 0.5: SNe within 30 deg of the mirrored
    major axis (17,109) have median r / R_e(semi-major) 0.79, of the minor
    axis (4,763) 0.71, ratio 0.90 where a truth ellipse would give ~q
    (~0.4). In elliptical radius 0.94 vs 1.84. Effective SN axis ratio
    ~q^0.8 (median-abs) / q^0.4 (rms). Measured vs truth size: moment
    major axis / truth 1.69 overall but 1.27 for resolved disks (R_e > 1",
    1.21 expected for exponential -> Gaussian), small hosts ~1.9 (PSF).
    Verdict: SNe sit on the mirrored ellipse at the truth radii along the
    major axis, but on a rounder footprint (about 2x the truth ellipse
    along the minor axis), not a smaller one; rendered disks are not much
    larger than truth. Connor's "SNe inside the visible galaxy" impression
    comes from the example selection: the 3" cap in baseline limits the
    largest-R_e examples (R_e 1.3-1.6") to offsets 1.3-1.8 R_e, and 72% of
    that pool is below 1.3 R_e; in sep10 the examples sit at 1.6-2.7 R_e.
  * *Operational.* lc remakes only on recipe/decision/input changes, never
    on a code change: changed scripts got explicit flags (`--tier dp02`,
    `--radius-page`, `--workers 16` on the ou24 candidate recipes) so the
    targets are remade. shared_interactive allows 2 jobs per user; a
    pending salloc that does not start within ~6 min is cancelled by lc
    ("cluster ... is ended"), and lc cluster names cannot contain
    underscores. A dropped ssh connection killed a `setsid nohup` driver on
    the login node; relaunch and lc resumes from the committed outputs.
- **lc on NERSC, operational.** `lc` hashes every declared input by content
  on every command with no cache (1 TB = 15 min), hence the objectTable index
  file; the 52 GB truth and 7 GB match directories cost ~1-2 min per command
  and are kept as real inputs. `lc compute launch` refused to start while
  `~/.lightcone` was group-writable (chmod 700). lc's result commits land in
  the NERSC clone only: `git push origin main git-annex` there before
  committing on the laptop, otherwise the branches diverge (a rebase of the
  DATALAD commits worked but left the annexed file "modified"; `git restore
  --staged --worktree results/` cleared it). Small results are brought to the
  laptop with scp + `git annex reinject <copy> results/<universe>/<file>` (the `--known` form reports "not known content" for lc-made files; the explicit two-argument form works, then `git annex smudge --update`).
- **DES tier.** 1623 Dovekie SNe Ia (IDSURVEY = 10) with HOSTGAL_RA/DEC from
  DES-SN5YR HEAD files; candidates from DES DR2 `des_dr2.main` via NOIRLab
  Data Lab TAP (use `q3c_radial_query(...) = 't'`; `phi = theta_j2000 + 90` in
  the +x West frame). DES computed DLR from A/B_IMAGE on the r+i+z detection
  image within 15 arcsec, d_DLR <= 4; reproducing the published host needs
  `dlr_radius: moments` and `search_radius: r15` in a DES-faithful universe.
  Prototype headline (Oct 2026, now rebuilt here): SBW recovered the DES host
  98.6%, DLR 98.4%; most residual failures were DR2 classing a faint host as a
  star, or different deblending.
- **DES tier (rebuilt, 29 Oct 2026, branch `des-tier`, laptop worktree
  `.claude/worktrees/des-tier-82bd51`; not merged; local compute only).** The
  published host is the DES-SN DLR assignment on deep SN-field coadds: every
  number is agreement with that reference, never recovery. Inputs:
  `data/des_sn5yr/` (Dovekie metadata + HEAD, sha256 in data/README.md) and
  `data/des_dr2_candidates.parquet` (655,978 unique DR2 objects within 2' of the
  1,623 SNe, fetched once 2026-10-29 by `src/des_fetch_dr2.py`, ~14 min over
  Data Lab TAP with 6 threads; not a recipe). `src/des_lib.py install()` swaps
  SExtractor light models into `dp02_lib` (shape = detection A, B, THETA_J2000
  for every band, written as `{b}_ixx/iyy/ixy` in arcsec^2 with zero PSF
  moments; exponential / sersic = n 1 / n 4 on the band's FLUX_RADIUS as
  determinant radius; moments = Gaussian on A, B; DLR moments = A, B; DLR
  half_light = FLUX_RADIUS ellipse; refExtendedness from extended_class_coadd
  >= 2). `src/des_run.py` runs dp02_associate / dp02_sweep (griz) / dp02_fig /
  dp2_fig unchanged and rewrites "true" -> "published" in every figure at save
  time (tick labels need their FuncFormatter wrapped), so no shared script was
  edited. Published host -> DR2 by `reference_match` in des_sample (99.14%
  matched at 1"). Random positions: 10 per SN in a 45-90" annulus. New universe
  `des_faithful` (dlr_radius moments, search_radius r15). All 17 des_* outputs
  current in baseline and des_faithful (~10 min total).
  * *Angle* (`des_angle_check`, 40 galaxies, r cutouts): theta_j2000 is East of
    North, <cos 2 delta> 0.926, median |delta| 3.7 deg (mirror -0.10, 55 deg).
  * *Agreement with the published host* (r, exponential, flb 0.01, d 4):
    baseline SBW argmax 97.78%, mean P_SBW(published) 0.973, DLR (half-light)
    97.10%, SBW = DLR 98.21%; des_faithful DLR (moments, 15") 97.60%, SBW
    unchanged. Prototype 98.6 / 98.4% does not reproduce exactly; the gap is the
    29 SNe (1.8%) whose published host is unmatched in DR2 (14; nearest DR2 object
    1.3-10" away, several with HOSTGAL_SNSEP 0, one 57" bright host) or a DR2
    star (15). Of the 1,594 whose published host survives the cuts SBW agrees
    99.56%, DLR 98.87% (faithful 99.37%).
  * *P-band matrix* (baseline): P > 0.5 / DLR = pub 1,575, / DLR != pub 8; P in
    [0.05, 0.5] / DLR != 6; P < 0.05 / DLR != 12; SBW hostless / DLR hostless 21,
    / DLR = pub 1. SBW and DLR fail together: the residuals are catalogue
    (absent / star), not method. HOSTGAL2 explains 44% of SBW-other and 73% of
    DLR-other picks; host-or-host2: SBW 98.77%, DLR 98.27%. ECE 0.013.
  * *Hostless*: random rejection at flb 0.01 75.8%, d 4 half-light 74.6%,
    d 4 moments 88.2%. Own J optima: flb 0.133 (J 0.858), d_DLR half-light 2.1
    (0.858), moments 2.7 (0.900). DP0.2 optima transfer at J -0.065 (flb
    0.0215) and -0.017 / 0.000 (d 2.7 half-light / moments); defaults lose
    0.114 (flb 0.01) and 0.125 / 0.033 (d 4). DES SNe sit close to their hosts
    (published median d_DLR 0.67), so the optimal flb is much higher than on DC2.
  * *Bands / profiles*: SBW host changes vs r < 0.6% in g, i, z; n 4 hostless
    0.06-0.18%, exponential 1.2-1.6%, moments 3.3%; argmax agreement 97.7-98.0%
    in every band and profile.
  * *Gallery*: sbw_ne_published cases are the published host classed as a DR2
    star or a knot / companion separated in the deep SN coadds but merged or
    absent in DR2. SIA footprints miss some tile edges (`des_gallery.py`
    re-searches 0.5 deg).
  * *Not done*: no `des_sample` decision (des_all_hosted needs ~12x the fetch);
    pos_3arcsec / frac_semimajor and ext_none universes not materialized for
    DES; a merge with dp2-legacy must add `dp2_candidate_catalogue` to
    `universes/des_faithful.yaml`.
- **DELIGHT tier.** 16,791 visually vetted hosts within 1 arcmin (so
  `search_radius: r60` for this tier); 14,111 ZTF/ALeRCE on Pan-STARRS DR1
  imagery plus 2,680 large nearby hosts (Asiago, HECATE). Candidate catalogue
  is DESI Legacy Surveys DR10 Tractor models (grz only, no i band; REX/SER
  types exist in DR10 though not in the Dey+19 overview). Check DR10 sky
  coverage against the DELIGHT positions before committing; Pan-STARRS stamps
  with the stamp form of SBW are the fallback (`delight_catalogue: ps1_stamps`).
  Labels are catalogue centroids that can be off by arcseconds, hence the
  `reference_match` decision. Keep the 2% of hosts with no SExtractor
  counterpart; DELIGHT kept them deliberately.
- **DP2 tier.** Samples live on NERSC under
  `/global/cfs/cdirs/lsst/groups/TD/SN/EDP2/SN_candidates/` (EDP2 transient
  WG; README, README_RC.md). Headline sample is `RC_sample.parquet` (690 with a
  good SALT3 fit, decision `dp2_sample: rc_salt`); alternatives are Mi Dai's
  `final_sample.parquet` (617), `parsnip_ddf.parquet` (7823 DDF ParSNIP) and
  `RC_lcfit_input_122k.parquet` (122k, minimal cuts). The laptop file
  `data/dp2_transients.csv` is exactly the parsnip_ddf id set (sha256 checked).
  `EDP2_Host_Associations/` there holds three runs on the 7823 set: J. Lee's
  DLR on Legacy DR10 Tractor shapes, Y. J. Kang's DELIGHT on cutouts, and
  Connor's earlier SBW csvs (8k and 122k); `dp2_cross_method_agreement`
  compares SBW with those. Candidates come from the DP2 Object table (both
  MultiProFit fits) over grouped TAP or the NERSC HATS `object_collection`.
- **DP2 clean-sample chain (28 Oct 2026, branch `dp2-tier`, NERSC clone
  `$PROJECT/host_association/paper-dp2`; not merged into main).** Inputs are
  index csvs: `dp2_transients` = `data/dp2_sn_candidates_index.csv` (the four
  EDP2 sample files; `src/dp2_sample.py` refuses to run if size/mtime moved)
  and `dp2_object` = `data/dp2_object_hats_index.csv` (the 8840 HATS
  partitions of `shared/rubin/DP2/HATS/object_collection/object_lc`; no
  lsdb: `src/dp2_candidates.py` reads whole partitions whose
  per_partition_statistics bounding box, padded by the radius, holds a
  position, so there is no margin cache). `src/dp2_lib.py` `install()` swaps
  MultiProFit light models into `dp02_lib` (exponential / free-n Sersic fits
  as the fitted Gaussian mixture, `SBW_CURVE`; DLR half_light = the
  exponential fit's half-light ellipse for every profile; moments as DP0.2)
  and `--tier dp2` in `dp02_associate/sweep/fig/gallery` uses it.
  `associate_one` gained `p_sbw_dlr` and a DLR-ellipse angle hook. New
  outputs: dp2_sample, dp2_candidates, dp2_association_sweep (all bands x
  profiles), dp2_random_positions (10 per transient in its own patch; DP2
  patches are ~10'), dp2_random_candidates/associations,
  dp2_hostless_indicator/metrics, dp2_p_top_host, dp2_disagreement_cases
  (+ ellipses), plus recipes for the declared ones. Materialized: baseline,
  physical (every DP2 file byte-identical to baseline except dp2_p_top_host,
  whose DP0.2 reference differs) and dlr_moments. Run cost: baseline ~40
  min on 16 CPUs, most of it lc hashing the DP0.2/OU2024 inputs; the
  recipes themselves take seconds to minutes.
  * *Sample (rc_salt, 690):* candidates 43,770 pairs within 30" (median 64
    per transient), ext05 removes 34% of them, none unusable.
  * *SBW vs DLR (r, exponential, flb 0.01, d 4 half-light):* same host
    85.65%, different 0.29% (2), SBW-only hostless 0.72% (5), DLR-only 0,
    both hostless 13.33% (92). P band of the DLR host: > 0.5 for 591/593
    hosted pairs. Near-identity is by construction: DLR half_light uses the
    same exponential ellipse SBW evaluates, so they differ only by flux
    weighting, and the 5 SBW-only hostless sit at d_DLR 3.7-4.0 (flb 0.01
    vs d 4 boundary). The informative comparison is `dlr_moments`: same
    83.91%, different 0.43%, SBW-only hostless 0.72%, DLR-only hostless
    1.59% (11), both 13.33%; there DLR-only hostless are mostly transients
    on a faint compact galaxy (flb 0.5-0.9) whose PSF-subtracted moment
    ellipse is tiny (d_DLR 5-18), plus a few where SBW's host is a large
    low-SB model.
  * *Hostless:* real kept / random rejected at flb 0.01 85.9 / 83.0%, flb
    0.0215 84.9 / 86.4%, flb 0.05 83.0 / 90.0%; d_DLR 2.7 82.3 / 90.5%, d 4
    86.7 / 81.3% (moments: 2.7 79.9 / 91.8%). DP2's own J optima are flb
    0.091 (J 0.733) and d_DLR 2.2 (0.734). The DP0.2 optima transfer at a
    cost of J 0.019 (flb 0.0215) and 0.006 (d 2.7); the defaults flb 0.01 / d
    4 lose 0.043 / 0.055. On DP2 SBW and DLR lie on the same ROC. 13-14% of
    the SALT sample is hostless for both (median nearest d_DLR 7.7; 25% at
    z > 0.6 vs 11-13% below). 99.7% of random positions still get a P > 0.5
    host, so hostlessness must come from flb, as on DP0.2.
  * *Top P:* SBW-hosted DP2 transients have median P_max 1.000, P > 0.9
    99.0% (DP0.2 97.8%), mean credible-set size at t 0.034 1.020 (DP0.2
    1.039): real-data SBW is at least as decisive as on DC2.
  * *Bands / profiles:* SBW host changes vs r in u 1.7%, other bands <=
    0.3%; DLR (half-light) 0.2% in u; hostless fraction u 40% (no u flux), g-y
    14-15%. Exponential vs Sersic same SBW host 100% (both hosted), vs
    moments 99.8%; the light models differ only in hostless calls (Sersic
    10.6%, exponential 14.1%, moments 19.3% SBW hostless).
  * *Gallery:* band flips are mostly within a deblended blend (siblings) or
    a large low-SB model winning in a faint band; sbw_ne_dlr are a compact
    object under the SN versus a larger neighbour's extended model.
    `src/dp02_gallery.py` now unwraps ObsCore patches straddling RA = 0
    (DP2 cutouts 422'd on one).
  * *Operational.* The laptop worktree shows annex-filter "modified" files
    that block `git pull --rebase`; push a wip ref, cherry-pick on NERSC,
    then `git reset --keep nersc/dp2-tier`. The DP0.2 inputs and
    dp02_associations content in the clone were reinjected from copies of
    the main clone's files (no `mainclone` remote). `pkill -f`/`pgrep -f`
    over ssh match their own command line.
    dp2_cross_method_agreement is still without a recipe (not started).
- **DP2 DELIGHT comparison (29 Oct 2026, branch `dp2-delight`, NERSC clone
  `$PROJECT/host_association/paper-dp2`; not merged).** Universe `dp2_ddf`
  (baseline + `dp2_sample: parsnip_ddf`, 7,823 DDF ParSNIP transients;
  557,091 candidate pairs, median 71 within 30"). Y. J. Kang's DELIGHT
  tables are the input `dp2_delight_yjkang` = `data/dp2_delight_yjkang_index.csv`
  (path, size, mtime, sha256; `src/ou24_index_files.py --sha256`); J. Lee's
  DLR and the fixed parsnip_ddf input are gone from the spec.
  `src/dp2_cross_method.py agreement|metrics|cases` matches the DELIGHT CNN
  position and the DELIGHT-SEx position to the transient's DP2 candidates
  (all, before cuts) under `reference_match` and reruns the association with
  the DELIGHT object as host to get its P_SBW and `true_host_status`; the
  SBW/DLR picks are asserted equal to dp2_associations. New outputs
  `dp2_cross_method_agreement` (csv), `dp2_cross_method_metrics` (json),
  `dp2_delight_cases`, `dp2_delight_gallery`, materialized only in dp2_ddf
  (other universes leave them unmade; the script would give the rc_salt
  overlap, 291 transients). Numbers (r, exponential, flb 0.01, d 4, 1"):
  * Hostless: SBW 52.9%, DLR 50.5%, DELIGHT 9.6% (no DP2 candidate within
    1" of its position), DELIGHT-SEx 1.6%. On rc_salt SBW is 14%: the DDF
    ParSNIP set is dominated by transients on compact sources.
  * Pairs (same / different / a-only hostless / b-only / both): SBW-DLR
    46.2 / 0.9 / 2.4 / 0 / 50.5; SBW-DELIGHT 27.4 / 14.3 / 48.7 / 5.4 / 4.2
    (same given both hosted 65.6%); DLR-DELIGHT 27.4 / 16.5 / 46.5 / 5.6 /
    4.0; DELIGHT-SEx 89.2% same.
  * The cause: DELIGHT's position lies a median 0.18" from the transient,
    and of the 7,072 DELIGHT positions that match a candidate, 4,762 (67%)
    are DP2 point sources (refExtendedness < 0.5) that ext05 removes before
    SBW/DLR; SBW calls 3,711 of those hostless and picks a neighbour
    (median 4.6" away) for 1,051. Where DELIGHT's object survives the cut,
    SBW picks the same object in 2,141 / 2,212 (96.8%); the 71 others have
    P_SBW(DELIGHT object) ~ 0. P bands of the DELIGHT host (n 7,072): > 0.5
    30.3%, t-0.5 0.1%, < t 15.7% (mostly the point sources, P = 0), SBW
    hostless 53.9%.
  * By image_survey: Legacy (5,100) SBW-DELIGHT same 28.3%, DELIGHT
    hostless 6.3%; PS1 (2,719) 25.7%, DELIGHT hostless 15.9% (shallower
    cutouts, more positions with no DP2 object); DES_DR2 has 4.
  * reference_match: 1" matches 90.4% of DELIGHT positions (Kang's own
    nearest-DP2 match, up to 2", 95.6%; the same object wherever both
    match), 3" 97.8% (SBW same 28.9%), frac_semimajor 16.0%.
  * Gallery: the transient sits on a compact, star-like DP2 source that
    DELIGHT and SEx both mark; SBW/DLR, having dropped it, go to a fainter
    extended neighbour or call hostless. A few DELIGHT/SEx splits are large
    nearby hosts where the CNN lands on the bulge and SEx on a knot. The
    open question for Connor is whether those point sources are the
    transient itself in the deep coadd, AGN/nuclear hosts, or stars
    (ext_none on DP2 would test it).
- **DP2 Legacy candidates (29 Oct 2026, branch `dp2-legacy`, NERSC clone
  `$PROJECT/host_association/paper-dp2` on that branch; not merged).** Decision
  `dp2_candidate_catalogue` (`dp2_object` default, `legacy_dr10`), input
  `legacy_dr10_sweeps` = `data/legacy_dr10_sweeps_index.csv` (14 DR10 south
  sweep tiles, 25 GB, within 0.3 deg of any parsnip_ddf transient;
  `src/legacy_index_sweeps.py`, sha256 checked against the survey's
  sha256sum), universe `dp2_ddf_legacy` (= dp2_ddf + legacy_dr10).
  `src/dp2_candidates_legacy.py` (called by `dp2_candidates.py --catalogue`)
  rewrites Tractor rows into the dp2_candidates schema so the whole DP2 chain
  runs unchanged: objectId = BRICKID * 1e6 + OBJID, refExtendedness 0 for
  TYPE PSF, DUP dropped, `exponential_*` = Tractor ellipse with n = 1,
  `sersic_*` = same ellipse with the Tractor n (REX/EXP 1, DEV 4, SER
  SERSIC, PSF 0.5), PSF sources a round Gaussian with R_e = PSFSIZE_R / 2
  (HWHM = Gaussian half-light radius), FLUX * 3631 = nJy for every flux,
  blendedness = FRACFLUX, no moments (moments profile / DLR moments are
  unusable on Legacy). Adding the decision to dp2_candidates,
  dp2_random_candidates and the two DP2 galleries makes every DP2 output in
  the other universes stale on this branch (identity change only; remake
  after merging). Galleries with legacy_dr10 draw the DP2 deep_coadd and the
  Legacy DR10 viewer cutout of each case side by side
  (`dp02_gallery.py --candidate-catalogue` / `--imaging rubin|legacy|both`).
  * *Tractor ellipse convention* (`dp2_legacy_angle_check`): SHAPE_R is the
    half-light *semi-major* axis (tractor `EllipseE.getRaDecBasis`), q =
    (1 - |e|)/(1 + |e|), and PA East of North = +0.5 atan2(SHAPE_E2,
    SHAPE_E1): on 27 large ECDFS galaxies the r-band cutout moments give
    <cos 2 delta> 0.77, median |delta| 3.9 deg for +theta, against 0.06 /
    40 deg for -theta. In the +x West frame phi = PA + 90 (as DES). The
    builder round-trips a, q, phi exactly through `sbw_host.sersic.ellipse_axes`.
  * *Kang's legacy_ls_id* is a string with the brick id at bit **20**, not 16
    (release << 40 | brickid << 20 | objid; bits 16-19 always zero): 99.7%
    then decode to a candidate at the identical position, and our 1" DELIGHT
    match is Kang's own Legacy counterpart for 7,097 / 7,117. Metrics recipe
    flag `--kang-ls-id-shift 20`.
  * *Coverage:* all four DDFs (ELAIS-S1, ECDFS, EDFS incl. Dec -48 to -51,
    COSMOS) are in DR10 south; 7,810 / 7,823 transients have a Legacy source
    with NOBS_R > 0 within 30" (median 47 vs 71 DP2 objects).
  * *Headline (r, exponential, flb 0.01, d 4 half-light, 1"):* hostless SBW
    72.0% (DP2 Object 52.9%), DLR 70.7% (50.5%), DELIGHT 9.0% (9.6%),
    DELIGHT-SEx 1.2% (1.6%). SBW-DELIGHT same / different / SBW-only
    hostless: 15.6 / 8.7 / 66.6% (DP2 27.4 / 14.3 / 48.7%); same given both
    hosted 64.1% (65.6%). DLR-DELIGHT same 15.1%, given both hosted 59.3%.
    SBW-DLR same given both hosted 96.5% (98.1%). P_SBW of the DELIGHT object
    (7,117 matched): > 0.5 17.2%, < t 9.4%, SBW hostless 73.2%.
  * *The DP2 point sources are not transient light.* Of the 4,748 DELIGHT
    objects that are DP2 point sources (1" match), the Legacy pre-transient
    catalogue has a PSF source at the same position for 4,704, a galaxy
    (REX/EXP/DEV/SER) for 17 and nothing for 27. Of those 4,596 distinct
    Legacy PSF objects (median r 19.6, IQR 19.0-20.2), 3,638 (79%) are Gaia
    stars (parallax S/N > 3 or proper motion S/N > 5) and 762 (17%) more
    have WISE W1 - W2 >= 0.8 (Vega; AGN colours); 146 other WISE-detected,
    50 undetected. The nearest Legacy source within 1" of the transient is
    PSF for 5,751 / 7,823 (74%). So the DDF ParSNIP set is mostly variable
    stars and AGN, and SBW / DLR calling them hostless after ext05 is right;
    DELIGHT, which always returns a host, marks the star / AGN itself. The
    hostless ROC agrees: real transients are barely distinguishable from
    random positions (Youden J 0.15 at flb 0.147 / d_DLR 2.0; at flb 0.01
    real kept 28.0%, random rejected 84.7%), unlike rc_salt (J 0.73).
  * *Hosted cross-tab (DP2 -> Legacy):* SBW hosted->hosted 1,769,
    hosted->hostless 1,915, hostless->hosted 424; DLR 1,839 / 2,033 / 450.
    Legacy loses hosts mostly through depth (fainter extended DP2 objects
    are absent or PSF in Legacy). Only 472 of the 4,748 point-source cases
    get an SBW host on Legacy, and in 14 it is the DELIGHT object.
  * *Gallery:* both cutouts show the same compact source under the
    transient before (Legacy) and during (DP2) the survey.
  * Numbers outside the outputs (point-source breakdown, Gaia / WISE,
    cross-tab) come from an ad-hoc script over the dp2_ddf and
    dp2_ddf_legacy results (session notes folder, `compare_legacy.py`); the
    ext05-free check of whether the remaining 958 non-Gaia PSF sources are
    AGN hosts was not done.
- **DP2 cleaned DDF sample (29 Oct 2026, branch `dp2-ddf-clean`, NERSC clone
  `$PROJECT/host_association/paper-dp2` on that branch; not merged).** New
  `dp2_sample` options computed in `src/dp2_sample.py` from `legacy_dr10_sweeps`
  (now an input of dp2_sample, so dp2_sample in every other universe is stale
  on this branch; only the new options read the sweeps): `parsnip_ddf_clean`
  drops transients with a Legacy DR10 source (not DUP) within 1" that is a
  Gaia star (REF_CAT G*, parallax S/N > 3 or PM S/N > 5) or AGN-coloured (PSF,
  W1 and W2 S/N > 3, W1 - W2 >= 0.8 Vega); `parsnip_ddf_galaxy` also drops any
  Legacy PSF source within 1". The sample carries gaia_star, wise_agn,
  legacy_psf_1arcsec, gaia_ref_no_astrometry, legacy_objectId, legacy_type,
  legacy_sep_arcsec, legacy_n_1arcsec, ddf_field, and the per-rule / per-field
  counts in `df.attrs["dp2_sample_accounting"]`. Universes `dp2_ddf_clean`,
  `dp2_ddf_clean_legacy`, `dp2_ddf_galaxy`; `dp2_ddf` gained its random /
  hostless outputs. All four materialized on one 4 h shared_interactive
  allocation in 1 h 26 min (~25 min per universe).
  * *Accounting* (sequential): 7,823 -> Gaia star -3,647 -> WISE AGN -1,575 =
    2,601 clean -> other PSF -546 = 2,055 galaxy. Per field (clean / galaxy):
    ELAIS-S1 2,238 -> 886 / 732, ECDFS 1,552 -> 496 / 401, EDFS 2,828 -> 726 /
    559, COSMOS 1,205 -> 493 / 363. Kept in clean though Gaia-referenced PSF
    without significant astrometry: 223. No Legacy source within 1": 902.
  * *Side by side* (dp2_ddf / clean / galaxy / ddf_legacy / clean_legacy; r,
    exponential, flb 0.01, d 4, 1"): hostless SBW 52.9 / 23.3 / 18.7 / 72.0 /
    38.7%, DLR 50.5 / 21.8 / 17.0 / 70.7 / 37.9, DELIGHT 9.6 / 16.7 / 17.5 /
    9.0 / 16.0, SEx 1.6 / 3.3 / 3.8 / 1.2 / 3.0. SBW-DELIGHT same 27.4 / 54.6 /
    56.8 / 15.6 / 46.5, different 14.3 / 10.8 / 12.2 / 8.7 / 6.2, SBW-only
    hostless 48.7 / 17.9 / 13.5 / 66.6 / 31.3; same given both hosted 65.6 /
    83.4 / 82.3 / 64.1 / 88.2 (DLR 62.4 / 81.9 / 80.6 / 59.3 / 84.1).
    P_SBW(DELIGHT object) > 0.5: 30.3 / 65.5 / 68.8 / 17.2 / 55.4%, < t 15.7 /
    12.6 / 14.5 / 9.4 / 6.9, SBW hostless 53.9 / 21.5 / 16.3 / 73.2 / 37.3.
    Hostless Youden J (SBW flb_best / DLR): 0.30 / 0.61 / 0.65 / 0.15 / 0.48
    (DLR within 0.005), optima flb 0.004 / 0.068 / 0.038 / 0.147 / 0.083; at
    flb 0.01 real kept 47 / 77 / 81 / 28 / 61%, random rejected ~82-85%. So
    the galaxy sample approaches rc_salt (J 0.73) but does not reach it.
  * *What is left* (clean, from an ad-hoc script `clean_diag.py` in the session
    notes folder): 634 of 2,167 matched DELIGHT objects are still DP2 point
    sources (galaxy: 436 / 1,695). sbw_hostless_delight_hosted 466, 426 of them
    a DP2 point source under the transient, with Legacy showing nothing (257)
    or an unclassified PSF (189); sbw_ne_delight 282, 208 with DELIGHT on a DP2
    point source and SBW on an extended neighbour whose 4 R_e ellipse covers
    the transient. delight_nocand 434: DELIGHT's position ~1.5-3" from the
    transient on no DP2 object while the transient itself sits on one (median
    0.16"); SBW hosts 295 of them. Gallery reading: the DELIGHT/SBW split is
    now mostly a faint compact source under the transient (point source in the
    DP2 deep coadd, absent or PSF in Legacy: compact hosts, nuclear transients
    or undetected AGN) versus a larger neighbour; where DELIGHT's object is
    extended SBW agrees almost always. Legacy candidates lose ~15 points of
    hosted transients to depth again.
- **Compute: NERSC Perlmutter via `lc`, set up 2026-10-26.** Lightcone runs
  `lc` *on the login node* and launches recipes into Slurm allocations from
  `~/.lightcone/compute.yaml` (local compute is refused on login nodes). On
  Perlmutter (`ssh cstone@perlmutter.nersc.gov`, alias `nersc`; project dir
  `$PROJECT` = `/global/cfs/cdirs/m1727/cstone`, account `m1727`):
  uv 0.13 and `lightcone-cli==0.5.0` (with git-annex) are installed in
  `~/.local/bin`; the uv cache must live on scratch (`UV_CACHE_DIR=$PSCRATCH/uv-cache`),
  so `~/.bashrc` needs `export PATH="$HOME/.local/bin:$PATH"` and that export.
  Two offers, both 16 CPUs / 60 GiB, one node: `shared-interactive`
  (`salloc`, qos `shared_interactive`, fast start, max 4 h) and `shared`
  (batch, qos `shared`, max 12 h). `lc compute launch --cpus 16 --memory 60
  --startup fast --wait` takes the first. This is by design: the analysis is
  catalogue queries and light numerics; raise the limit only if a recipe needs
  it (deblending explorations later may). The project clone on NERSC is
  `$PROJECT/host_association/paper`, with the bare hub
  `$PROJECT/host_association/paper.git` as git remote `nersc` on both sides
  (`git push nersc main git-annex` from the laptop, pull on NERSC). Annexed
  data does not yet transfer through the hub (the laptop cannot reach the hub's
  git-annex-shell over ssh): copy the file with `scp` to `$PSCRATCH` and run
  `git annex reinject --known <file>` in the NERSC clone, as was done for
  `data/dp2_transients.csv`. Perlmutter has no landlock, so `lc` runs recipes
  there unsandboxed (`sandbox: none`) but still records every run. The laptop keeps its own
  `local` offer for probes over TAP.
- **Rejected / parked.** OU2024 with our own Sersic fits on its coadds was
  rejected on the wrong belief that only one preview coadd patch survives;
  revisit with the OU2024 tier (above). DP2 DiaObject SN-like
  selection as the DP2 sample (parked: the collaborator list was preferred).
  Pareto-style "best method" tables (out of scope by framing).
- **Prototype code** for DES, DP0.2 and OU2024 from 22 Oct 2026 no longer
  exists under `host_association/`; the user chose to rebuild cleanly here as
  Lightcone outputs. `../hostassoc-repo/examples/` still has the DP2 demo and
  disagreement-search scripts.
- **Provenance.** 63 prior insights from 8 papers (Sullivan+06, Gupta+16,
  Qu+24, Gagliano+21 GHOST, Forster+22 DELIGHT, DESC DC2, OpenUniverse2024,
  Dey+19 Legacy Surveys), all quotes verified against cached PDFs on
  2026-10-26; IOP DOIs would not download, so the arXiv DOIs are cited.
- **Merging branch results into the main NERSC clone (29 Oct 2026).** Merging
  dp2-delight, dp2-legacy, dp2-ddf-clean and des-tier into main brought their
  `results/` into `$PROJECT/host_association/paper` as annex pointers (the hub
  stores no annex content). `lc status` reports such outputs *current* (the
  manifest hashes match) while a recipe that reads them fails on an empty
  pointer ("Parquet magic bytes not found", "lacks columns"). Before
  materializing on a clone that did not make an output, give it the content:
  small files by `git annex reinject --known <copy>` from another clone or the
  laptop, large ones by key with
  `git -C <making-clone> annex contentlocation <key>` and the same reinject.
  Index csvs in `data/` need the same treatment. DES data (`data/des_sn5yr/`,
  `data/des_dr2_candidates.parquet`) and all `des_*` results live on the laptop
  only; on NERSC they show stale but are not meant to be made there. After the
  `dp2_sample` / `dp2_candidates` recipe changes, remade tables hashed
  identically, so only `dp2_random_candidates` and the galleries were actually
  rebuilt; the rest came back current by content hash.
- **Paper configuration = spec defaults (30 Oct 2026; spec edited, nothing
  materialized yet).** `dlr_radius` default is now `moments` (the literature
  DLR, independent of SBW's ellipse); `dp02_sample_cuts` stays `physical`.
  `universes/baseline.yaml` is exactly the spec defaults (physical cut +
  moment DLR); every other universe is baseline plus its one change, except
  the two Legacy-candidate DP2 universes, which keep `half_light` because
  Tractor carries no moments. The appendix universes keep their historical
  ids because a universe id is part of every rendered recipe (`results/<id>/`
  appears in the manifest's recipe string), so renaming a results directory
  makes all of its outputs stale: `physical` = baseline with `half_light`
  DLR (its decisions are unchanged), `dlr_moments` = baseline with `clean3`
  (unchanged). The DELIGHT training-set tier (inputs `delight_hosts`,
  `delight_candidate_objects`, outputs `delight_*`, decision
  `delight_catalogue`, `search_radius: r60`) is removed: the Förster+22 set
  is not public. `summary_table` now takes `ou24_recovery` too. With the
  new defaults `baseline`, `sep10`, `truth_nearest`, `ext_none`,
  `des_faithful` and the DP2 universes are stale by decision change and
  need a remake; `physical` and `dlr_moments` are untouched. The plan and
  the audit behind this are in `OUTLINE.md`.
- **Calibration is candidate-level (30 Oct 2026; code and spec changed, not
  yet materialized).** `calibration_curve(df, cut)` in `dp02_lib.py` bins
  every candidate with P > P_TOP_MIN (the `p_top` lists, float32) of every SN
  in `calibration_subset` (SBW keeps it as hosted, flb_best >= cut, and the
  true/published host is a surviving candidate) and compares the fraction of
  candidates per bin that are the true host with the bin's mean P. The
  calibration error is sum_bins |n_true - sum P| / N_SNe. The old curve (top
  host only) is gone; `dp02_fig.py calibration`, `dp02_sweep.py` (band /
  profile panels, profile_table) and `recovery.json`'s
  `expected_calibration_error` all use the new one, and the calibration and
  sweep recipes carry `--hostless-cut {decisions.sbw_hostless_cut}`. Probed
  on DES baseline: error 0.012 per SN on 1,589 hosted SNe. DP0.2 / OU2024
  figures will be made with the baseline remake.
- **Sample ladder in the tables; miss accounting rewritten (31 Oct 2026;
  code and spec changed, remake pending).** The simulated samples
  (`dp02_sample`, `ou24_sample`) are rung 2 of the ladder in `OUTLINE.md`
  section 2: hosted, host row present, *detectable* under the new decision
  `transient_peak_mag` (default `r24`: true peak r < 24 from the new input
  `dp02_sn_peak`, the DC2 SN truth summary, with t0 inside the `dp02_visits`
  window 59582-61406; OU2024 uses SNANA `PEAKMAG_r`), within the offset cut.
  That keeps ~10 % of hosted DC2 SNe (~41k) and ~4 % of OU2024 transients
  (~21k); the counts at every step are in the parquet metadata key
  `sample_ladder` and are copied into `*_miss_accounting.json["ladder"]`
  (the faint fraction is a number for the text, never drawn). SNe whose true
  host has no catalogue object are *kept* with `host_objectId = -1`
  (`true_host_status = absent`), so the association tables are the whole
  rung 2 and every consumer filters its rung: `dp02_lib.rung(df, 4)` (true
  host kept by the object-table cuts and with a DLR ellipse) is what every
  truth figure and metric uses (matrix, ROC, calibration, recovery-vs-X,
  P(true), credible sets, offsets, band/profile sweeps, in DES too: 1,594 of
  1,623); `recovery.json` has `ladder` and `by_rung` (rungs 2, 3, 4). The
  miss accounting (`dp02_miss.py`) draws rung 2 and, inside cat1 and cat2,
  a bar-in-bar of the pipeline's outcome for SNe whose host it could not
  find (hostless call vs wrong host assigned; `outcome` in the json);
  panel (c) is gone. The old `n_unmatched` bookkeeping is gone with it.
  DES outputs were already stale from the calibration change, so the rung
  filter reaches them in the same remake.
