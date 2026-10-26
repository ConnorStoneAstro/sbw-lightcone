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
- **Simulation tier is DP0.2, not OpenUniverse2024.** The public OU2024 release
  has no multiband Object table (only skyCatalog truth, single-visit src, one
  preview coadd patch); DP0.2 has CModel bulge+disk fits plus SN truth with
  host ids over the RSP TAP (`dp02_dc2_catalogs`). On DP0.2 "exponential" maps
  to the CModel disk-only component and "sersic" to bulge+disk; the stamp form
  needs deblender model stamps from the butler. OU2024 truth is used only for
  the cosmology propagation. Both simulations place SNe by sampling the host
  light profile (DC2: Sersic n = 1 or 4 component chosen by stellar mass;
  OU2024: SNANA HOSTLIB surface-brightness weight), so SBW is favoured by
  construction there; say so in the paper. OU2024 hosts are also ~1 dex too
  low in stellar mass (known issue).
- **DP0.2 quirks** (from the 22 Oct 2026 prototype): ~5% of SN truth rows point
  to an empty host row; ~9% sit > 3 arcsec from the declared host; all z > 1
  and 10% of z < 1 SNe are deliberately hostless (control sample, not
  failures). Keep objectIds int64 end to end; sync TAP joins time out, use
  chunked `id IN (...)` queries or async.
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
- **DP2 tier.** `data/dp2_transients.csv` holds 7823 positions supplied by a
  collaborator (18-digit DiaObject-like keys); the selection behind the list is
  undocumented and must be obtained from the provider before publication.
  Candidates come from the DP2 Object table (both MultiProFit fits) over
  grouped TAP or the NERSC HATS `object_collection` via LSDB.
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
  (`git push nersc main git-annex` from the laptop, pull on NERSC; annexed
  data moves with `git annex copy --to nersc`). The laptop keeps its own
  `local` offer for probes over TAP.
- **Rejected / parked.** OU2024 with our own Sersic fits on its coadds
  (rejected: only one preview coadd patch is public). DP2 DiaObject SN-like
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
