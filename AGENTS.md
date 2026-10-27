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
  user decided OU2024 is folded back in as a second simulated tier; that
  work is being prepared in a separate session, and the spec does not carry
  it yet. The OU2024 prior insights are still in the spec as evidence.
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
  * *Radius fallback* (peer request): 10 SBW misses in the cutout prototype
    were hosts with NaN `r_bdReD` (p = 0), so a galaxy ~20" away took P = 1.
    `light_model(..., radius_fallback=True)` now falls back to the other CModel
    radius, then the moment half-light radius. The association and sweep
    recipes gained `--radius-fallback`, so the *recipe identity changed* and
    everything downstream was remade in all five universes. Baseline: SBW
    97.42 -> 97.52%, both-wrong-same 4,221 -> 4,177, DLR unchanged.
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
    matchers. Class sizes: truth_matchers_disagree 75, sbw_wrong 7,591,
    dlr_wrong 26,846, sbw_ne_dlr 22,740, both_wrong_same 4,737 (no flb cut).
    `src/dp02_gallery.py` is `../dp02_cutouts_proto/gallery.py` copied
    unchanged (mtime 2026-10-27 11:20 EDT, sha256 e6c7f5a3...880b2; that
    directory is not a git repo). `dp02_disagreement_gallery` has a recipe
    but is *not materialized*. The laptop sandbox reads the declared token
    `rsp_token` (`../.rsp_image_tap_token`) and reaches data.lsst.cloud
    (checked with `lc run`; a probe run made 59/60 cutouts in 90 s, 5 pages).
    But `lc materialize` on the laptop wants to rebuild the whole chain,
    because the NERSC-only inputs cannot be checked there. The options are to
    put the token at `$PROJECT/host_association/.rsp_image_tap_token` on NERSC
    (same relative path, mode 600; compute nodes have internet), or to
    materialize from a laptop clone that has those inputs. That is the user's
    call, since it copies a credential.
  * *Operational.* The NERSC clone must `git fetch origin && git annex merge`
    before `git push origin main git-annex`, or the git-annex push is
    rejected. A remade output whose bytes are identical leaves its downstream
    only "behind", not stale. `lc compute launch ... --time 4h` gives the full
    shared_interactive window. The full five-universe rematerialization with
    random chains took 1 h 50 min on 16 CPUs. `$PROJECT` paths inside
    `n.sh`-style double-quoted ssh commands expand on the laptop, so quote
    remote globs in single quotes.
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
