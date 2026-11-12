# Changelog

## 2026-11-12b — FACTS land-water storage dropped past 2100 (Marcus)

- `python/figures/figlib.py`: `FACTS_CLIMATE_DRIVEN` loses "ssp-lws" at 2150 and 2300. FACTS's land-water module reads
  no climate (SSP population, held at its 2100 value), so its post-2100 High-minus-Very-Low difference was a population
  difference; it had been kept as an exception. Same change as the research repo (SLR-RFF-BRICK CHANGELOG 2026-11-12d).
  The comparison table loses 14 rows (546 remain) and the responsiveness table 2 (78). `check_figures_identity.py`
  passes 16/16 against the re-rendered research figures.

## 2026-11-12 — Fig. 5: short panel titles (Marcus)

- `python/figures/plot_vv_gsic_ladrillo.py`: the panel (b) and (c) titles drop their explanatory clauses ("the 3
  reservoirs equilibrate, so declining T slows the melt"; "falls toward zero as the reservoirs equilibrate"), as the
  research script does from the same day (SLR-RFF-BRICK CHANGELOG 2026-11-12). `check_figures_identity.py` passes 16/16
  against the re-rendered research figure.
- Later the same day, panel (a) drops "(FaIR ensemble mean)" (the caption says it); identity gate 16/16 again.

## 2026-11-11f — Post-processing: the four production chains give the shipped posterior, byte for byte

The second of the four calibration-port pieces (SLR-RFF-BRICK `notes/handoff_2026-11-11c_table4_seeds_calibration_scope.md`
§5). The target builders and the Code availability text are still to come.
- New `src/postprocess.jl` (module `Ladrillo.Postprocess`) and `scripts/postprocess_calibration.jl [--chains=] [--outdir=]
  [--accept-slr] [--force]`. They port the research post-processing of 2026-09-20 (`run_l27_postprocess.sh`): its
  noise-mode check (second-half median sd_gis < 0.10 cm in every chain), `diag_slr_convergence_by_chain_ladrillo.jl`
  (R-hat and ESS of projected SSP2-4.5 total sea level at 2100 and 2150, 400 draws per chain) and
  `postprocess_mcmc_ext.jl --tag=L27 --accept-slr` (R-hat and ESS of every parameter with the ESS maxlag capped at
  200,000; acceptance; every 400th pooled draw; the pooled covariance). The research ran these as three programs that
  passed the sea-level table through a file guarded by a freshness check; here they run in one process.
- The sea-level test runs on the projection model as it stood when L27 was accepted: land water `:central`, Greenland tap
  off, (lambda, T_crit) attached in one block of 400 per chain (seed 20260920, so the four chains get the same rows) and
  gamma at its paleo median. These are named constants in `Postprocess`. The shipped projections use `:observed` land
  water and the `:single` paleo assignment.
- Outputs: the posterior in the format of `data/posterior/`, `slr_convergence_L27.csv`, `parameter_convergence_L27.csv`
  (every parameter; the research log printed only the failures), `adapted_cov_L27.csv` with the parameter names as
  header (the research file had x1..x50), and a provenance file with the seed and where it is applied.
- Not ported: the research script's six-median comparison against the 2018 BRICK baseline subsample (a file the package
  does not ship), and the check of each chain's log banner (replaced by checking the chain header against the parameter
  list; the research SLOWP/FAST names are read as SLOWG/FASTG).
- New direct dependency MCMCDiagnosticTools 0.3.18, the version both Manifests already held. `Project.toml` was edited by
  hand and `Pkg.resolve()` changed only the Manifest's project hash and the package's dependency list.
- **Gate** (`test/calibration/gate_postprocess.jl <SLR-RFF-BRICK repo>`). PASS on the reference runtime. It checks the
  four production chains' sha256, runs the script on them with `--accept-slr` (25 min, 21 of them the parameter ESS;
  3.3 GB peak) and requires (1) the posterior byte-identical to `data/posterior/ladrillo_L27_posterior.csv`, (2)
  `slr_convergence_L27.csv` byte-identical to the research file (R-hat 1.0007869032040153 at 2100, 1.0016940768436582 at
  2150), (3) the covariance equal to the research file at all 2,500 entries, and (4) 28 log lines equal to the research
  log: the noise check, the per-chain projected quantiles, the eight Antarctic parameters that fail (R-hat up to 1.313,
  ESS 19.5-243.8) and the acceptance lines.
- `tools/build_test_reference.py`: new kind `log_lines` (the research log's lines matching the patterns the gate
  repeats). A research source that is not under version control is now accepted and recorded as `untracked` in
  `SOURCES.csv`, where its sha256 pins it. Three new references: `research_slr_convergence_L27.csv`,
  `research_adapted_cov_L27.csv`, `research_postprocess_L27_log_lines.txt`.
- `test/runtests.jl`: 118/118 pass, unchanged. The gate stays outside the suite because it needs the 8 GB of production
  chains.
- **Mutation test** (`test/calibration/mutation_test_postprocess.sh`, run on a frozen copy of this commit; result in
  `mutation_result_postprocess.txt`). Four mutants are killed: P1 subsample offset (check 1), P2 certificate thinning
  (checks 2 and 4), P3 covariance jitter 1e-10 to 1e-9 (check 3, 50 of 2,500 entries) and P4 ESS maxlag 200,000 to
  100,000 (check 4). The no-power control P5 (noise-mode threshold 0.10 to 0.20 cm) survives, as it should: every
  chain's second-half median sd_gis (0.025 cm) is far below both thresholds.

## 2026-11-11e — The calibration sampler; it reproduces the production L27 chain

**Marcus 11-11 (session 2):** port the calibration drivers and target builders before the release. This entry is the
first of the four pieces scoped in SLR-RFF-BRICK `notes/handoff_2026-11-11c_table4_seeds_calibration_scope.md` §5: the
calibrator and its gate. Pooling the chains and the target builders are still to come.
- New `src/calibration.jl` (module `Ladrillo.Calibration`): the L27 log-posterior, the parameter table, the start rows,
  the proposal-covariance seed and `run_chain`. It ports the research `julia/calibrate_mcmc_ext.jl` with the flags of
  `run_mcmc_L27.sh` (43d84aa), reduced to that configuration: no flags, no other arms. The model is this package's
  `LadrilloModel` on 1850-2026 (land water `:central`, Greenland tap off). Every likelihood and prior term keeps the
  research arithmetic and order of accumulation.
- New `scripts/calibrate.jl <n_iter> <seed>`: one chain, its adapted covariance and a provenance file recording the seed
  and where it is applied (`Random.seed!` on the global RNG immediately before `RAM_sample`).
- New `data/calibration/` (imported by `tools/import_l27_data.py`, label-only renames): the four over-dispersed start
  rows, the L26 adapted covariance that seeds the proposal, and the Greenland amplification prior table.
- New direct dependencies Distributions, LinearAlgebra and RobustAdaptiveMetropolisSampler. `Pkg.add` first upgraded
  Distributions 0.25.126 -> 0.25.134; reverted and pinned, so the Manifest now equals the research `julia_v2`
  environment on every shared package.
- **Gate** (`test/calibration/gate_chain_prefix.jl`, new references from the four production chains via
  `tools/build_test_reference.py`; their sha256 in `test/reference/SOURCES.csv`). PASS: the log-posterior at all four
  start rows equals production (792.9208399273703 for seed 2026), and the first 20,000 rows of the seed-2026 chain equal
  production in the log_post of every row and in all 50 parameters at 44 reference rows. The window holds 4,603 accepted
  moves; the research gate's 300 rows held none (the first move is row 460). 20,000 iterations take 103 s here
  (production ran at ~140/s).
- **Mutation test** (`test/calibration/mutation_test.sh`, result in `test/calibration/mutation_result.txt`). Three
  mutants leave the log-posterior at all four start rows unchanged and are killed by the chain alone: the sampler's
  target acceptance 0.234 -> 0.235 (row 460), the AR(1) rho bound 0.99 -> 0.98 (row 460; the start rows' largest rho is
  0.976) and the active-basin share targets written as 0.799 / 0.816 (row 796, a 4e-13 difference in log_post). The
  research gate's 300 rows would miss all three. The GlaMBIE partition sd 0.05 -> 0.0501 is killed at the start rows. The
  no-power control (the geometry block's diagonal proposal scale, which the L26 covariance seed overwrites) survives.
- `test/runtests.jl` gains the start-row check and a check of the parameter table against the research calibrator's own
  `--dump-priors` output (`data/constants/ladrillo_priors_L27.csv`): all 50 names, components, means, sds and bounds equal.
- **The production chain held gamma at 2.038, not at its paleo median.** The first port, which held lambda, T_crit and
  gamma at their paleo medians as the research flags `--cut-fastdyn` / `--fix-gamma` describe, missed the four
  production start values by -0.015 to +0.034. A state dump against the research calibrator (every glacier, Greenland,
  thermal-expansion and land-water series identical; Antarctica not) found the cause: the research calibrator sets the
  paleo medians at :1605-1607 and then `update_brick_nu3!(m, medoid, ...)` at :1627 loads the central row over all three.
  Holding the central row's values (lambda 0.01365, T_crit -15.14, gamma 2.038) reproduces all four starts exactly;
  holding the paleo lambda and T_crit with the central gamma does too, so lambda and T_crit are inert in the hindcast
  and gamma is the only one that matters. `posterior()` attaches gamma = 2.834 (the paleo median) to every product.
  Paired on 500 draws (rows 1:20:10000, ensemble-mean forcing), 2.834 against 2.038: hindcast AIS differs by at most
  0.015 cm (AIS RMSE 0.1276 vs 0.1274 cm); projected AIS by at most 0.13 cm through 2100; in 2300 the AIS median is
  lower by 1.4 cm (SSP2-4.5), 3.2 cm (SSP5-8.5) and 2.6 cm (vvH), and the SSP5-8.5 total median by 4.1 cm.
- Would a calibration at gamma = 2.834 differ? Importance weights on all 10,000 draws, w = exp(logpost at 2.834 minus
  logpost at 2.038), every other term identical: the log-posterior difference spans -0.27 to +0.42 (p05-p95 -0.07 to
  +0.10), the effective sample size is 9,972 (99.7 %) and the largest weight is 1.5 x uniform. Reweighted p05/p50/p95 of
  AIS and the total at 2100 and 2300 (SSP2-4.5, SSP5-8.5, ensemble-mean forcing, gamma 2.834) move by at most 0.54 cm
  (SSP2-4.5 2300 total p95, 407.87 -> 407.33), most by under 0.1 cm. So the shipped products, which run gamma 2.834,
  match what a calibration at 2.834 would give to about half a centimetre. The port holds what production ran
  (`Calibration.HELD`); how the release describes this is Marcus's decision (README *Known issues*).

## 2026-11-11d — Every file a figure reads records the seeds behind it

**Marcus 11-11 (session 2):** yes to recording seed provenance in the figure inputs, so the paper can say where the
seeds are recorded.
- **Audit.** The projection, hindcast and BRICK 2.0 files already recorded their seeds: the paleo-draw assignment
  (`PALEO_SEED` 20260920), the draw-to-config pairing (2026), BRICK 2.0's `Random.seed!(2026)` before `get_model`, and
  the posterior predictive's noise (`NOISE_SEED` 2026). Every `rand`/`randperm`/`randn` call in `src/` and `scripts/` is
  one of these. Two gaps:
  - `outputs/tables/vv_model_comparison.csv`, which Figs. 2, 3 and 6 read, recorded only the pairing seed.
  - No provenance said whether the FaIR forcing carries internal variability. The SSP files and the calibration forcing
    ssp245harm (Fig. 1's hindcast) keep the calib 1.6.0 posterior's `stochastic_run=True` with `use_seed=True` (841
    configs, 841 distinct seeds in its `seed` column); the van Vuuren files are deterministic.
- `Ladrillo.fair_seed_note(scenario)` (`src/forcing.jl`) states which, and the provenance of `project_joint.jl`,
  `brick20_joint.jl`, `project_fixed.jl`, `posterior_predictive.jl`, `posterior_predictive_brick20.jl` and
  `ic_hindcast_residuals.jl` now carries it.
- `vv_model_comparison.py` copies every provenance clause that names a seed from the 14 draws files into its own
  provenance, and stops if a model's clauses differ between scenarios or a required seed (paleo and pairing for
  Ladrillo; `Random.seed!` before `get_model` and pairing for BRICK 2.0) is missing. Mutation-tested on scratch copies:
  a paleo seed changed in one file, BRICK 2.0's seed clause dropped from one file, and the Ladrillo pairing clause
  dropped each stop the script; the unmutated copy passes.
- **Two runs are byte-identical.** `posterior_predictive.jl`, `posterior_predictive_brick20.jl`, and `project_joint.jl`
  and `brick20_joint.jl` on vvLN, each run twice: all 11 output files byte-identical between the runs.
- **The local products were regenerated** (both hindcasts and the seven van Vuuren scenarios for both models, 53
  files) in a staging directory, checked equal to the previous products on every non-provenance cell (0 differences),
  and swapped in. The SSP projection products were not regenerated and keep their 10-08 provenance (no figure reads
  them). `check_figures_identity.py`: PASS, every PNG still byte-identical to the paper's render.
  `check_tables_identity.py`: PASS. `test/runtests.jl`: 97/97, exact.
- Not changed: the research repository's figure inputs (`vv_model_comparison_L27.csv` and the
  `scope_slr_fairunc_{draws,paths}` files have no provenance column). Recording seeds there means re-running the research
  projections and changing column sets that this package's identity gates read.

## 2026-11-11c — Table 4's script; Table A1 shows its log₁₀ rows in ln, as the research table does

**Marcus 11-11 (session 2):** the Code availability text claims the scripts for every table, so Table 4's scorecard
moves into the package now; the calibration drivers and target builders move before the release.
- New `python/tables/ladrillo_vs_brick20_scorecard.py`, a port of the research
  `python/scope_ladrillo_vs_brick20_scorecard.py`: mean bias, RMSE, largest error and 90 % band coverage of Ladrillo and
  BRICK 2.0 against the hindcast observations, by component and in five windows. It reads the two posterior-predictive
  band files in `outputs/hindcast/`; `--ladrillo`, `--brick20` and `--targets` let the test suite point it at the
  research references. Its provenance records that it uses no random numbers.
- Its inputs were checked first: both band files equal the research files on every non-provenance column, and the
  targets file is the research `recalib_targets_ext.csv` (md5 070f74ab).
- `check_tables_identity.py` gains a Table 4 row: 250 of 250 numbers bit-identical to the research
  `outputs/scope_ladrillo_vs_brick20_scorecard_L27.csv`. New references `research_postpred_brick20_components_timeseries.csv`
  and `research_scope_ladrillo_vs_brick20_scorecard_L27.csv` (SLR-RFF-BRICK 88ae6b4); `test/runtests.jl` runs the
  scorecard on them, so the suite checks its arithmetic from a fresh clone.
- **Table A1 had fallen behind the research table.** On 11-10 (research b8d7fdf, Marcus "Option 1") the research
  builder began showing the four parameters the calibrator samples as log₁₀ (the three glacier κ and `gis_s_high`) in
  natural logs, with footnote ᶜ, so that Table A1 uses one log base. The package port still printed log₁₀; the gate
  against the research repository failed on five lines of the markdown. Ported: the markdown scales prior μ, σ, bounds
  and the posterior quantiles by ln 10, and the [LN-DISPLAY] check recomputes each row as quantiles of ln(10^x). The
  CSV keeps log₁₀ values and units, as the research CSV does. References for both Table A1 files re-copied from the
  research repository.
- `mutation_test_tables.py` adds seven mutants (Table A1's prior left unscaled; ln 10 written as 2.30; Table 4's window,
  band column, one observation and one BRICK 2.0 median moved by 1e-12, and its re-reference window). The two that must
  stop the script (the ln 10 value and the re-reference window) do; the other five fail the gate. It also fixes the
  harness's handling of a script that stops (it returned two values where three were unpacked). Result in
  `python/tables/mutation_result.txt`: 14 of 14 as expected.
- Verified: `test/runtests.jl` 97/97, exact on the pinned runtime; `check_tables_identity.py` against the research
  repository PASS on all five outputs.

## 2026-11-11b — Version 1.0.0; the internal builds become v0.1–v0.4; the v0.1 options renamed

**Marcus 11-11 (GMD guidance, Decision 1):** the paper reports Ladrillo v1.0. The internal tags are renamed and the
release is tagged v1.0.0. **Ruled the same day:** the options that reproduce the first build get the flag `--v0.1` and
version-free option names; push the commits, hold the tag changes on GitHub until the release.

| internal build | old tag | new tag | commit |
|---|---|---|---|
| v0.1 (first build) | v1.0.0 | v0.1.0 | a2bb3b0 |
| v0.2 | v1.1.0 | v0.2.0 | 78b03c6 |
| v0.3 | v1.2.0 | v0.3.0 | c526331 |
| v0.4 | v1.3.0 | v0.4.0 | 01462a8 |

Entries below this one keep the version names they were written with.
- `Project.toml` and `CITATION.cff` are 1.0.0; `date-released` and the `v1.0.0` tag wait for the release commit. The
  README's version text and *Versions* section use the new numbers and add a v1.0 entry.
- Renamed: `--v1.0` → `--v0.1` (the three projection scripts), `posterior(paleo=:v1)` → `paleo=:per_selection`,
  `lws_anchor=:v1_step` → `lws_anchor=:obs_at_first_year`; output suffixes follow (`_paleoper_selection`,
  `_lwsobs_at_first_year`). `--v1.0` now stops with an error rather than silently running the defaults.
- Comments that called the first builds v1.0–v1.3 now say v0.1–v0.4; text about the model or the L27 posterior, which
  every version shares, no longer carries a version.
- Identity gate: `cases.jl` uses the package's names; `dump_old.jl` translates them to the research kernel's own
  (`RESEARCH_NAME`), whose environment variables and defaults (`v1_step`, `v1`) are unchanged.
- The tags were renamed in this clone only. GitHub keeps v1.0.0–v1.3.0 until the release.
- Verified: `test/runtests.jl` 96/96, exact, reporting Ladrillo 1.0.0. The old option names now fail with an error
  (`paleo=:v1`, `lws_anchor=:v1_step`, `--v1.0`). The research side's translation was checked by re-dumping K3 and
  K3_v11 through it: byte-identical to the 10-08 reference dumps.
- **Identity gate 1 PASS, bit-identical on all eight cases, after one stale reference was replaced.** The first run
  passed K1–K3 and failed K4 by up to 19.9 cm. The K4 reference dumps were written at 06:19–06:21 on 10-08, before
  the vvH cube was rebuilt on the published emissions that afternoon (the v0.3 change), and gate 1 was not re-dumped
  for K4 then; v0.3 and v0.4 were certified by gates 2–4. The research kernel re-run on today's cube (sha256 equal to
  `data/forcing/`) is bit-identical to the package on both K4 cases. The stale dumps are in
  `test/identity/out/quarantine/20261111_k4_stale_vv_forcing/` (local; `out/` is not tracked) with a README.

## 2026-11-11 — The paper's figure scripts (Figs. 1–6) and the FACTS / MAGICC-SLR comparison data

**Marcus 11-11 (GMD guidance, Decision 3):** the Code availability text says the archive holds the scripts that
produce every figure and table, so the figure scripts move into the package before the release.
- `python/figures/`: ports of the six research builders (`plot_hindcast_components.py`, `vv_model_comparison.py`,
  `plot_model_comparison_components.py`, `plot_future_components.py`, `plot_vv_gsic_ladrillo.py` from
  `plot_vv_gsic_wr_vs_ladrillo.py --ladrillo-only`, `plot_vv_responsiveness.py`) and `figlib.py` (cut down from
  `ladrillo_figs.py`). Drawing code copied unchanged; only the I/O differs. Paper mode only, one posterior, outputs to
  `outputs/figures/` and `outputs/tables/`.
- Inputs are the package's own products. Before porting, every file the figures read was compared with its research
  counterpart: the joint and BRICK 2.0 paths, cells and draws for all ten scenarios, both posterior-predictive files,
  the calibration targets and the mean-GMST drivers are equal on every non-provenance column. The package's gates
  files have no CONTROL row, so `figlib.gate_ladrillo` requires PAIRING, SPLICE-IDENTITY and SUM instead.
- `data/comparison/`: FACTS (n = 200) and MAGICC-SLR outputs for the van Vuuren scenarios, MAGICC's GMST, and the
  FACTS module classes, copied unchanged from SLR-RFF-BRICK, with four `data/MANIFEST.csv` rows and a README.
- Dropped from the research scripts: the opt-in overlays of Fig. 1 (MAGICC history, IGCC GMSL, IMBIE, the >2000 m band)
  and their input files; the SSP sets and memo mode; the `_width`, `_gmst` and `_pathdep` side tables (no figure
  reads them); the Wigley–Raper glacier file behind Fig. 5, whose panel (a) GMST now comes from
  `data/forcing/fair_mean_gmst_*.csv` (equal on every drawn row).
- Identity: `python/figures/check_figures_identity.py` passes 16/16 against the research renders in
  `figures/paper/`: all seven PNGs byte-identical, all captions identical except Fig. 5's driver-commit clause
  (declared), the comparison table (560 × 15) and the responsiveness table (80 × 11) equal. Mutation-tested
  (`python/figures/mutation_result.txt`): a one-pixel change, a caption edit, a 1e-12 cell change, a stale declared
  difference and a missing PNG each fail.
- `requirements.txt` pins matplotlib 3.10.8, the version the paper renders were made with.

## 2026-10-15c — Table 5's .md: computed reading sentences and a versioned label

- `python/ic_ladrillo_vs_brick20.py` (6674807, c964ccb):
  - The steric-row and ρ-at-bound sentences are computed, as in SLR-RFF-BRICK (its CHANGELOG 2026-10-15h). "A tie by
    construction" was false on the published BRICK 2.0 posterior: Δln L (steric) is +7.2 / −0.8 / −2.3 at
    ρ ≤ 0.99 / 0.95 / 0.90.
  - The label was a literal "Ladrillo v1.0" through v1.3.0. It now reads the version from `Project.toml`, as the 10-09
    banner fix does.
- Re-run at all three ρ bounds. Every CSV and perdraw CSV is byte-identical to its pre-run copy, and identity gate 5
  passes. Only the .md files change. `outputs/tables/` is not tracked, so the v1.3.0 tag is unaffected.

## 2026-10-09 — Readiness review for sharing with Tony; v1.2.0 re-tested from a GitHub clone; stale labels fixed

**Asked (Marcus 10-09):** is the GitHub repo ready to share with Tony?
- **Fresh clone of `v1.2.0` FROM GITHUB** (`gh repo clone`, default depot; the empty-depot `Pkg.instantiate()` of
  10-08i still holds, because the Manifest changed only in its self-entry since): `test/runtests.jl` **90 pass + 2
  broken, exit 0, every comparison EXACT**, 7 min 07 s at load ~35. The 2 broken checks are the BRICK 2.0 posterior's
  MANIFEST row and comparison; that file is not shipped. The README quick start runs (10,000 draws, 451 years; the
  components sum to the total: 51.83 cm at 2100 on draw 1, SSP2-4.5).
- **Tony's 10-07 request had two parts.** "Obtain Ladrillo from GitHub and run an example" is met. "Rerun the full
  workflow" is not yet: the calibration code (target prep, the MCMC, the postprocess) and the figure scripts are not
  ported, and neither are the static FACTS / MAGICC-SLR comparison data the figures read.
- **Fixed (labels and comments only; no code path changes):**
  - `scripts/posterior_predictive.jl` and `scripts/ic_hindcast_residuals.jl` printed a hard-coded "v1.0 (L27)"
    console banner at v1.2.0. It now comes from `pkgversion(Ladrillo)`, like the files' provenance columns, which were
    already right.
  - `src/record_crossing.jl` said "the paper (v1.1) does not condition", which contradicts the v1.2 release. It now
    says the paper's joint arm conditions each draw on its own climate, in `scripts/project_joint.jl`.
  - README layout: `data/posterior/` holds the L27 posterior, unchanged since v1.0 (it had said "the v1.0 posterior");
    `data/brick20/` is fetched, not shipped. The *BRICK 2.0 posterior* section says where to put the file if you
    already have it.
  - `test/identity/fix_effects.jl`: dropped a reference to a private configuration file.
- `test/runtests.jl` after the fixes, in the working checkout (BRICK 2.0 posterior present): **92 / 92, EXACT**, 5 min 45 s;
  the posterior-predictive banner prints `1.2.0 (posterior L27)`. The `ic_hindcast_residuals.jl` banner is not in the
  suite; it was checked on its own (it prints `Ladrillo 1.2.0 (posterior L27) | 2000 of 10000 draws | ...`).
- **Still open before PUBLIC** (10-08h/j list, re-checked): the BRICK 2.0 posterior URL (Tony); the AR6 / Rennels
  file's origin and licence (2 of 43 forcing files); the calibration and figure ports; a data-licence and attribution
  section (the README states only the code licence); `~/climate-env` as the tests' Python default; this CHANGELOG is a
  development log full of internal rulings; a Zenodo archive of the release for GMD.

## 2026-10-15b — v1.3.0 released: the version the GMD paper reports

**Marcus 10-15:** "update the version # for Ladrillo." v1.3.0 = v1.2's model and numbers + the BRICK 2.0 arm on the
published Zenodo posterior (2026-10-15 below).
- `Project.toml` and the Manifest's self-entry: 1.2.0 → 1.3.0.
- `CITATION.cff`: version 1.3.0, date 2026-10-15; its note names v1.3.0 as the paper's tag.
- README: the intro and *Versions* state v1.3 and what it changed.
- No code change since ac1ccf0, which certified the products bit for bit.
- Tagged `v1.3.0`; `main` and the tag pushed. The repo stays PRIVATE.
- **Clean clone of the tag from GitHub** (`git clone --branch v1.3.0`, then `tools/fetch_brick20_posterior.sh`, which downloaded
  the posterior from Zenodo and verified its sha256): **92 / 92, nothing skipped**, BRICK 2.0 612 / 612 bit-identical. The
  v1.2.0 clone test skipped 2, because the posterior could not be fetched.

## 2026-10-15 — The BRICK 2.0 arm reads the published MimiBRICK v2.0.0 joint posterior from Zenodo

**Ruled (Marcus 2026-10-15):** "move the BRICK 2.0 arm to the Zenodo joint file." The arm now reads the 35 BRICK
columns of `parameters_subsample_sneasybrick.csv` from Zenodo 10.5281/zenodo.20592337 (CC-BY-4.0, sha256 `d74577f6…`).
Through v1.2.0 it read Tony Wong's May-22 subsample (sha256 `78247db8…`), which had no public URL.
- `tools/fetch_brick20_posterior.sh` downloads it and checks the sha256. The blocker in 10-08h and i is gone: a fresh
  clone can now fetch the BRICK 2.0 posterior.
- `src/brick20.jl`:
  - `BRICK20_FILE` names the new file, and `BRICK20_LABEL` is the provenance label the three scripts now use.
  - The file stores `antarctic_precip0` in log space, so `apply_brick20!` no longer logs it.
  - `brick20_posterior()` refuses a file whose precip is positive in every row (a linear-units file).
- `data/MANIFEST.csv` swaps the row; `tools/import_l27_data.py` points at the new source.
- `python/ic_ladrillo_vs_brick20.py` (Table 5) parses the header as CSV and does not charge the 16 SNEASY columns.
  K stays 27 physical + 8 noise, and the assertion is unchanged.
- `test/reference/research_brick20_joint_draws_ssp245.csv` was rebuilt from the re-run research output
  (`tools/build_test_reference.py --only`, SLR-RFF-BRICK 96da010). 510 of its 612 values moved; the rest are land water.

**Before the reference was rebuilt,** the suite gave 91 pass / 1 fail, and the fail was the BRICK 2.0 comparison
(510 values outside tolerance). After: 92 / 92, BRICK 2.0 612 / 612 bit-identical.

**Identity gates 2-5 PASS, bit-identical against the research re-run** (SLR-RFF-BRICK CHANGELOG 2026-10-15d):
products for all 10 scenarios (Ladrillo and BRICK 2.0), panel, posterior predictive, hindcast (`postpred_brick20`, both
residual files, obs sigma), and Table 5 at rho 0.99 / 0.95 / 0.9. Ladrillo's residuals and obs sigma differ from the
v1.2.0 outputs only in the provenance version string ("1.2.0-DEV" → "1.2.0").

**What moved:** BRICK 2.0 numbers only. SSP2-4.5 fixed arm, BRICK 2.0 total 70.66 / 317.72 → 67.75 / 304.18 cm at
2100 / 2300. Table 5 ar1_prof ΔAIC 88.9 → 95.5. The full list of paper numbers is in SLR-RFF-BRICK
`outputs/brick20_swap_paper_number_diff_20261015.csv`.

**Quarantine:** the v1.2.0 BRICK 2.0 products are in `outputs/quarantine/20261015_brick20_may22_posterior/` (local;
`outputs/` is not tracked).

**Version:** `Project.toml` still says 1.2.0, but `main` now differs from tag `v1.2.0` in the BRICK 2.0 posterior. The
next tag is Marcus's call.

## 2026-10-08l — v1.2.0 released: the version the GMD paper reports

**Marcus 10-08: the conditioned numbers are "Ladrillo v1.2".** v1.1 cannot reproduce them: the same inputs give the
unconditioned joint arm and the harmonized-tail vv forcing.
- `Project.toml` and the Manifest's self-entry: `1.2.0-DEV` → `1.2.0`.
- `CITATION.cff`: version 1.2.0. Its note now names v1.2.0 as the paper's tag and v1.1.0 as the unconditioned
  predecessor.
- README: the intro states v1.2 and what each version changed. *Versions* drops "in development". *Known issues* no
  longer says the label is undecided.
- No code change since 2a1a34b, which certified the conditioned products bit for bit (10-08k).
- `test/runtests.jl` at 1.2.0: **92 / 92, EXACT** (reference runtime; it prints `Ladrillo 1.2.0`). 4 min 18 s at load ~48: a CCX job
  shared the machine.
- Tagged `v1.2.0`; `main` and the tag pushed. The repo stays PRIVATE: making it public is Marcus's step.

## 2026-10-08k — Joint arm conditioned on the record by default; the vv re-import and the conditioning certified bit-for-bit

Marcus 10-08: "go with b and any pending fixes that can be done simultaneously". His rulings:
- the vv forcing moves to the published emissions;
- the joint arms take the record-conditioned draws ("own-config, drop");
- MAGICC-SLR is re-run;
- the pulse arc waits.

The research side is SLR-RFF-BRICK CHANGELOG 2026-10-08k.

- **`scripts/project_joint.jl`:**
  - Drops every joint draw whose Antarctic fast dynamics fires by `record_end()` (2025) on the GMST it actually runs
    with (`fastdyn_onset_year` on its own spliced config).
  - Statistics and the draws file use the kept draws; the [SPLICE-IDENTITY] and [SUM] gates still see all draws.
  - The gates file lists every rejected draw with its onset year: SSP draws 229, 610, 1260; vv draws 187 (2024) and
    610 (2022).
  - `--no-record-conditioning` gives the v1.1 joint arm, suffixed `_uncond`; `--v1.0` does not condition.
- **Certified** (`test/identity/run_products.sh`, 24 runs, then `compare_all.sh` against SLR-RFF-BRICK f7942da+):
  - Ladrillo joint arms: 3 SSP × 71,946 values and 7 vv × 71,964 values, all identical. The row counts are 2,000 minus
    3 (or 2) dropped draws × 18, so the two codebases dropped the SAME draws.
  - BRICK 2.0: 10 × 36,000 identical.
  - Panel, posterior predictive, BRICK 2.0 hindcast and the IC residuals: identical.
- `test/runtests.jl`: 92 / 92, exact.
- README *Versions* and *Known issues* now describe v1.2 as it is: new vv forcing and the conditioned joint arm. The
  paper's version label for these numbers is Marcus's call.

## 2026-10-08j — The vv forcing rebuilds from public inputs (a names-only skeleton replaces the IIASA template); an absent input no longer blocks the build

Follows 696a603 (the vv files re-imported on the published Zenodo emissions). Closes items 1, 2, 5 and 6 of the
10-08h pre-public list for `tools/`, and item 4 except the first fetch.
- **Skeleton instead of the IIASA template.** On the `zenodo` variant every value from 2023.5 on comes from Zenodo
  20713982 v1.1.1, and the cmip7harm builder reads only the 2023.5 anchor and later columns. The template contributed
  only names, units and its column grid.
  - New `tools/forcing/vv_template_skeleton.csv`: 357 rows (7 markers × 51 species, template row order), columns
    `slug, scenario, region, variable, unit`, no numbers.
  - `build_vanvuuren_extended_splice.py --skeleton` (only with `--variants zenodo`) reproduces the template grid,
    including its gap (2021.5 and 2022.5 absent), and fills every column from Zenodo × `EXT_TO_SMITH_SCALE`. Its
    history columns are therefore Zenodo's own history, not Smith 2024; nothing downstream reads them. A species
    absent from Zenodo is an error. `--spliced-dir` is unchanged for the other variants.
  - `build_forcing.py`: `VV_TAIL_VARIANT = "zenodo"`; step 05 (`build_vanvuuren_spliced.py`) removed from the plan;
    step 06 uses the skeleton and depends on nothing. The script stays for the template-based variants.
- **Inputs 15 → 13.** `vv_iiasa_prerelease` and `smith2024_history` removed from `forcing_inputs.INPUTS`,
  `inputs.sha256` and the fetch notes (grep: step 05 was their only reader). The AR6 entry stays; step 03 needs it.
- **An absent input no longer blocks the build.** Each step declares its inputs. Argv inputs are recorded as the
  plan is assembled. The calib steps' raw files come from the setup scripts' own constants (`NEED` etc.). Step 04
  declares the three RCMIP files that `fill_from_rcmip` reads.
  - ABSENT input: the steps that read it and everything downstream are dropped, and the build prints which steps
    and which of the 43 files are skipped, then builds the rest.
  - PRESENT but failing sha256 (RCMIP: or fair's md5): the build aborts before any step.
  - `--install` refuses unless no step was dropped and all 43 exist.
  - `check_forcing_identity.py --allow-missing` gates only what exists and lists the absent files; without the flag
    nothing changes.
  - A plan self-check: the declared outputs must equal the 43 expected files, and every INPUTS key must have a reader.
  - `--intermediates` now compares `spliced_ext_zenodo` on its names and its 2023.5+ columns.
- **Provenance text.** README rewritten for the chain without step 5, the 13 inputs and the skip behaviour. Licences
  were re-read from the APIs (metadata only): Zenodo 18828694 CC-BY-4.0 (attribution added to the README,
  `forcing_inputs.py` and the `fetch_inputs.py` printout), Zenodo 20713982 MIT, chrisroadmap/fair-calibrate MIT;
  RCMIP not checked. The fair-calibrate URL is pinned to `3f72545` (the last commit touching the file, from
  `gh api .../commits?path=`); nothing was downloaded, so its bytes are first checked by `fetch_inputs.py --verify`.
  The root README's "two inputs are not public" line is corrected.
- **Internal references trimmed** in `tools/` comments, docstrings and printed strings (no code change): memory and
  skill names, two research branch names, a session id, handoff paths, the NCC-revision lines, "Marcus's machine",
  personal paths in usage text, and stale `scripts/` paths. The BRICK 2.0 posterior lines no longer describe a
  pending arrangement. Two generated calibration READMEs (intermediates, not `data/forcing/`) now name
  `tools/forcing/`.

**Gates (2026-10-08, load average ~35, `--jobs 2`):**

| gate | result |
|---|---|
| (a) full rebuild, `--fairtofredi-root`, all inputs (179 s) | **43 / 43 sha256-identical** to `data/forcing/` (identity gate and a separate `cmp`); 7 / 7 `emissions_v160_cmip7harm_vv<M>.csv` (and their `.variant`) byte-identical to the research copies; `--intermediates` 0 FAIL (23 CSVs byte-identical, 7 splices equal from 2023.5 on) |
| (b) `--inputs-dir` copies of every input except the AR6 file (147 s) | steps 03 and 08 dropped and named; **41 built, 41 / 41 sha256-identical** (`--allow-missing`: 41 PASS, 2 ABSENT); not built: `fair_mean_{gmst,ohc}_ssp245harm.csv`; the default mode still FAILs the 2 missing; `--install` refuses |
| (c) mutation | one digit changed in a copy of the Zenodo vv input: abort before any step (exit 1, no logs), both with all inputs and with AR6 absent. One digit in a built `fair_cube_ohc_vvML_raw.csv`: FAIL (max \|diff\| 5.0e-4), with or without `--allow-missing`. +1 ulp in a 2150.5 splice cell: `--intermediates` FAIL (2.8e-17). Gate self-test PASS |
| (d) `git grep` for the five private strings under `tools/` | no match |

**What a reviewer can rebuild from public inputs alone: 41 of the 43 forcing files**, all 28 vv files included.
The two ssp245harm files need the non-public AR6 file.
- Still open: the BRICK 2.0 posterior has no public URL. The AR6 file has no recorded origin or licence. The
  fair-calibrate URL's bytes are unverified until first fetched. The package's vv certification (gate 2, per the
  696a603 note) was not re-run here.

## 2026-10-08i — The suite FAILED from a fresh GitHub clone (90 / 1 fail); the BRICK 2.0 posterior's MANIFEST row now skips when absent

**Found** by the fresh-clone test FROM GITHUB (README's literal `git clone`, empty `JULIA_DEPOT_PATH`):
- `Pkg.instantiate()`: 3 min 18 s, every dependency downloaded and precompiled from scratch, MimiBRICK from its
  `v2.0.0` tag. ⚠ `Pkg.status()` prints MimiBRICK as
  `v1.2.0-dev`: that is the version field in MimiBRICK v2.0.0's own Project.toml, not a wrong pin.
- Quick start: runs, 10,000 draws, 451 years.
- `test/runtests.jl` and `Pkg.test()`: **90 passed, 1 FAILED, 1 broken**. The failure was the data-MANIFEST check on
  `brick20/parameters_subsample_brick.csv`, which ruling 6 says is fetched, not shipped, so EVERY fresh clone lacks it.
  The 92 / 92 of 10-08e/g ran only in Marcus's checkout, which has the file; so did the "Pkg.test() also passes" check.
  The README's "from a fresh clone; nothing outside this repository is needed" was false.

**Fix (Marcus: skip-if-absent):** `FETCHED_DATA` in `test/runtests.jl` names that row, derived from
`Ladrillo.BRICK20_FILE`. It is checked when the file is present and recorded as `@test_skip` when it is absent, as
the BRICK 2.0 testset already is. The README's Validation section says so.

**Checked on three frozen copies of the fixed tree:**

| state | result |
|---|---|
| file absent | 90 pass, 2 broken, exit 0 |
| file present, correct | 92 / 92, exit 0 |
| file present, one digit changed (byte 498) | **MANIFEST check FAILS**, exit 1 |

⚠ In the third state the BRICK 2.0 comparison still PASSES. The changed digit is in row 1, which IS sampled (rows
1:10:10000), but in column `sd_glaciers`, a likelihood noise parameter that never enters a model run. (First guess, that
the row was unsampled, was checked and was wrong.) For the noise columns, and for the 9,000 unsampled rows, the checksum
is the only guard.

**Re-run from GitHub at ab55bd2** (a fresh clone, using the depot the empty-depot run built): `test/runtests.jl` gives
90 pass + 2 broken, exit 0, every comparison EXACT; `Pkg.test()` gives 90 pass + 2 broken, exit 0, in portable mode,
max relative difference 1.7e-15. Tables A1/A2 used `~/climate-env`'s Python; without NumPy and pandas those checks
skip as well.

## 2026-10-08h — Private GitHub repo `msarofim/Ladrillo.jl` created; `main` and tags `v1.0.0`, `v1.1.0` pushed

Per the 10-07 ruling: Claude creates the repo PRIVATE; Marcus makes it public and adds Tony. Name confirmed by Marcus
(`Ladrillo.jl`, the Julia convention).
- **Pre-push checks** (the repo will later be public, so they were run as a publication review):
  - working tree clean at 08a9b99; refs are `main` and the two tags only; no stash;
  - no path ever committed is absent from HEAD, so the 10-07 purge of the BRICK 2.0 posterior left nothing behind;
  - no `/Users/`, scratch, HPC, email or token strings in any tracked file or in any added line of history;
  - `test/runtests.jl`: 92 / 92, EXACT, 2 min 06 s.
- **Commit identity.** All 15 commits are authored `Marcus Sarofim <MarcusMarcus@MacBook-Pro-2.local>`, because the
  repo had no `user.email`. GitHub cannot link that address to the account.
  - Marcus: keep the history as it is (a rewrite would change every SHA, including both tags, and both repos'
    CHANGELOGs, the notes and memory cite them).
  - New commits use the repo-local identity `Marcus C. Sarofim <msarofim@users.noreply.github.com>`.
- **Remote:** tags resolve to v1.0.0 = a2bb3b0 and v1.1.0 = 78b03c6, as certified.
- `CITATION.cff` gains `repository-code`; the README's Install section gets the clone URL. Both 404 for anyone
  else until the repo is public.

**Before the repo is made PUBLIC** (a delegated read of all 22 `tools/forcing/` files plus the three other `tools/`
scripts; the rebuildability claim was re-checked here). None of these blocks the private push:
1. ⛔ **The 28 shipped `data/forcing/*vv*` files derive from the gated IIASA `scenariomip-2026-prerelease`**, whose
   redistribution terms are not established. They are DATA in the repo, not just a missing input. Options: clear it
   with IIASA, wait for the final public ScenarioMIP release, or rebuild on Zenodo 20713982 (changes numbers; the
   `zenodo` variant still reads the step-5 template, which itself needs the IIASA file).
2. **The forcing README overstates what a reviewer can rebuild.** The 10-08b handoff said a reviewer could rebuild
   the SSP and calibration forcing. Through `build_forcing.py` they can rebuild NONE of the 43 files:
   - `build_forcing.py` calls `fi.verify()` over all 15 inputs whatever `--stages` says, and stops on any failure;
   - `build_emissions_v160_ssp245.py` requires `--ar6-file` unconditionally, so `ssp245harm` (the calibration driver)
     needs the Rennels file, not only its post-2100 tail.
   The 13 SSP files that need only public inputs can be rebuilt by running steps 1, 2, 4, 9 and 11 by hand. Fix:
   verify only the inputs the requested stages read, or say so in the README.
3. The EDF / Rennels AR6 file has no recorded public origin or licence.
4. Licences now checked through the APIs: Zenodo 18828694 is CC-BY-4.0 (attribution required); Zenodo 20713982 and
   both chrisroadmap repos are MIT; RCMIP is unchecked. The chrisroadmap `fair-calibrate` URL points at the moving
   `main` branch; pin it to a commit. Nothing has been downloaded and compared against `inputs.sha256` yet.
5. Internal references to trim: memory and skill names, the branch names `heat-ed-morbidity` and
   `claude/calib160-migration`, session ids, private handoff paths, the "NCC revision / submitted version" lines in
   `build_vanvuuren_extended_splice.py`, "on Marcus's machine" / `~/Documents/...` paths in the forcing README,
   `build_forcing.py`, `check_forcing_identity.py`, `forcing_inputs.py` and `import_l27_data.py`, and the
   collaborator-status lines in `fetch_brick20_posterior.sh` and `import_l27_data.py`.
6. Stale `scripts/...` paths in five forcing builders (the files now live in `tools/forcing/`). Two are written into a
   generated README; none reaches a forcing CSV, so correcting them is safe.
7. `test/runtests.jl` defaults to `~/climate-env` for Python, and `requirements-forcing.txt` names it too.
8. Unchanged from 10-08: the BRICK 2.0 posterior has no public URL.

## 2026-10-08g — Tables A1 and A2 ported (`python/tables/`); the test suite gains them (92 checks)

Ported by a delegated agent and verified here.
- **Data:**
  - `data/constants/ladrillo_priors_L27.csv`, the calibrator's prior dump with the glacier labels renamed.
  - `data/constants/paleo_geo_prior_ton.csv`, a byte copy.
  - Both come through `tools/import_l27_data.py --only`, with a new `labels` transform under the existing
    reversibility check. Two MANIFEST rows are appended; the 55 existing rows are byte-identical.
- **Code:** the three research scripts (`diag_ais_block_pca.py`, `ladrillo_table_a2.py`,
  `ladrillo_prior_posterior_table.py`), numerical code copied.
  - The research script's own `_FAST → _FASTG` label map is DROPPED: the names arrive renamed, and mapping them again
    would print `_FASTGG`. The gate checks for `FASTGG`.
- **Gate** (`check_tables_identity.py`, its rename regex written independently of the import tool):
  - the PCA CSV: 266 / 266 numbers bit-identical;
  - Table A2 `.md` and Table A1 `.md`: text identical, apart from the provenance comment;
  - Table A1 CSV: 473 / 473 numbers bit-identical, and 12 text cells equal after the rename.
- **References:** test/reference gains all four research outputs (SOURCES rows appended, sha-checked copies).
  - **Tried and fixed:** the agent's version copied only Table A2's `.md` and Table A1's CSV. Table A2 prints 3
    significant figures, so the self-contained suite could not see a 1-ulp input change. Adding the PCA CSV closes
    that.
  - Checked: a 1-ulp change to one paleo-prior sd (`ais_slope`) now FAILS the reference-only gate (235 PCA values),
    while the Table A2 markdown stays identical.
  - ⚠ Its largest move is 5.5e-13 relative, inside the portable rtol, so only the EXACT run sees it.
- **The agent's mutation run** (`python/tables/mutation_result.txt`, against the research repo):
  - a 1e-12 posterior perturbation, a 1-ulp prior sd, PROP_MIN either way and a doubled rename all FAIL;
  - a 1-ulp posterior perturbation correctly survives, because pandas' default float parser never delivers it to the
    calculation.
- **`test/runtests.jl`: 92 / 92** in 2 min, exact.

**For Marcus:**
1. The research table scripts read the posterior with pandas' DEFAULT float parser, which is not exact, and the ports
   keep that, so the tables reproduce. A re-implementation with an exact parser will differ in the last bits.
   Invisible at the tables' precision. Same family as the Table 5 comparator fix of 10-07e.
2. Table A1's heading still says "SLOWG = SLOWP and FASTG = FAST in the code". It is kept word for word so the table
   matches the paper. It no longer describes this package's code; the wording is Marcus's.

## 2026-10-08f — The FaIR forcing builders ported (`tools/forcing/`): all 43 shipped forcing files rebuild byte-for-byte

**Item 3** (Marcus: "3, 5, and 6 should be done first"). Ported from the private FaIRtoFrEDI @ 8984a13 by a delegated
agent and verified here.
- **Contents:** the calibration setup, the emissions builders, the mean / cube / spliced / van Vuuren cube builders,
  and two van Vuuren upstream splice builders the chain needed. Only the I/O plumbing changed; the numerical code is
  copied, and the diff against the originals is plumbing and docstrings.
- **Entry and gates:**
  - Entry point: `build_forcing.py --work-dir W {--fairtofredi-root R | --inputs-dir I} [--jobs N]`. It writes
    `data/forcing/` only with `--install`.
  - Every input is checked against `inputs.sha256`; nothing is ever downloaded mid-build. `fetch_inputs.py` holds the
    download code.
  - `check_forcing_identity.py` is the gate.
- **Verified independently in this session:**
  - a fresh rebuild (55 s, `--jobs 3`) → 43 of 43 PASS in the gate;
  - a separate `shasum -a 256` comparison of every rebuilt file with `data/forcing/` → 0 differ.
  - The agent reports four rebuilds (`--jobs 1/3`, both input layouts) all byte-identical, every intermediate
    byte-identical to the research copies, and the gate self-test passing.
- **Mutation:**
  - the calibration forcing with the obvious `calibration_v160 --marker none` → GMST FAIL (−0.030 K at 2024, max
    0.037 K), matching memory `l24_forcing_rebuild_command`;
  - vvM with `--stochastic on` → FAIL.
- **⛔ Two inputs are NOT public:**
  - The IIASA `scenariomip-2026-prerelease` harmonized + infilled emissions (gated login, redistribution not
    established) feed ALL 28 van Vuuren files, 2023–2100.
    - Checked as a substitute: the public Zenodo 20713982 v1.1.1 `emissions_1750-2500.csv` is NOT one. Only 14,196 of
      the IIASA file's 30,030 rows match by species name (CO2 is split differently), only 18% of those are equal, and
      minor gases differ by up to 4× at 2100.
  - The EDF / Rennels `AR6_emissions_ssp245_1750_2300.csv` (no licence recorded) shapes ssp245harm after 2100 only,
    outside the calibration window.
  - The build stops at its input check without them. Marcus decides.
- `calibration_v160_prod/` (research-only, no LICENCE) is no longer an input: the setup steps regenerate it
  byte-identically from Zenodo 18828694.

## 2026-10-08e — A self-contained test suite; Julia compat 1.12; CITATION 1.1.0; the BRICK 2.0 URL question

**`test/runtests.jl`** (Marcus: "3, 5, and 6 should be done first"). It runs from a fresh clone, with no research repo.
- **References:** `tools/build_test_reference.py <SLR-RFF-BRICK>` copies values the RESEARCH drivers wrote into
  `test/reference/` (264 KB). Lines are copied as text, so no float is re-printed. It refuses dirty sources.
  `SOURCES.csv` records each file's source, commit (SLR e2d8827) and sha256.
- **Tests (83 checks):**
  - every data file against `MANIFEST.csv`;
  - 33 joint-arm draws × both arms × SSP2-4.5 and SSP5-8.5 (2,376 values). The set includes the three draws that cross
    the fast-dynamics threshold inside the record, and checks the pairing itself;
  - the full posterior predictive (bands, bias, coverage);
  - the component-sum and splice identities;
  - the v1.0 options (paleo moves only λ / T_crit; the land-water step moves glaciers / GIS / TE not at all, and land
    water by exactly the step);
  - the record-crossing check, as a subprocess;
  - 612 BRICK 2.0 values, skipped with a notice if its posterior is absent.
- **Exactness:**
  - On the REFERENCE RUNTIME (macOS arm64, Julia 1.12.6, default bounds checking) every comparison is `==`: all 2,988
    reference values are bit-identical.
  - Elsewhere the comparison is to rtol 1e-12 (~451 steps × eps × 10).
  - **Tried and fixed:** the first `Pkg.test()` FAILED 4 checks that the direct run passed. Test: a direct run with
    `--check-bounds=yes` reproduced exactly the same counts (241 / 183 / 21 / 4). So it is Pkg.test's forced bounds
    checking re-associating float operations: ≤ 1.7e-15 relative, coverage counts unchanged, only `mean_bias` moves
    by 1–2 ulp. The reference runtime now requires default bounds checking. `Pkg.test()` runs in portable mode and
    passes.
  - **Tried and fixed:** `Printf` was missing from the test targets.
- **Mutation test** (scratch copies; `src/` untouched):
  - 1 ulp on `AIS_TANT0` (the Antarctic temperature map, a live path) → 4 testsets fail.
  - 1 ulp on BRICK 2.0's `te_s₀` → the BRICK testset fails (13 values).
  - 1 ulp on `GIS_V0_M` SURVIVES. It enters only clamps that never bind (`_gis_eq`'s [0, v0] clip, the whole-sheet
    headroom), the same no-power control as `mutation_test.sh`'s.
  - Portable mode cannot see a 1-ulp mutation by construction; only the exact run can.

**Item 6:**
- `[compat] julia = "1.12"`, was "1.11" though the Manifest is 1.12-format and the README says 1.12.
- `Pkg.resolve()` prunes NOTHING from the Manifest: all 306 entries are MimiBRICK v2.0.0's own dependency closure
  (Plots, XLSX, NetCDF, …), which cannot be trimmed without changing MimiBRICK. Its stale self-entry version
  (1.0.0 → 1.2.0-DEV) is refreshed; no other entry changes.
- `CITATION.cff` → version 1.1.0 (the paper's tag), date-released 2026-10-08, with a note that `main` may run ahead.
- `[extras]` / `[targets]` for the tests: Printf, SHA, Test.

**The BRICK 2.0 posterior URL: not found.**
- MimiBRICK v2.0.0's README points to Zenodo 10.5281/zenodo.20592337 (2026-06-08, CC-BY-4.0). Its
  `parameters_subsample_brick.csv` (7,011,352 B, md5 `6487623e…`, sha256 `fecabef3…`) is a DIFFERENT calibration from
  ours (6,621,725 B, sha256 `78247db8…`, Tony's post-PR#93 subsample of 2026-05-22):
  - 0 of 350,000 cells match;
  - `antarctic_precip0` is stored in log space;
  - the median `thermal_alpha` is 0.054 against 0.158.
- The three earlier versions of the record (17966962, 6626335, 6461560) do not match either.
- Like-for-like on this package's FaIR forcing (1,000 draws each):

  | | TE hindcast RMSE 1900–2025 | SSP2-4.5 TE 2100 | other components |
  |---|---|---|---|
  | ours | 0.52 cm | 19.4 cm | — |
  | Zenodo | 1.97 cm | 6.6 cm | within ~15% of ours |

- So it is not a drop-in substitute. My guess, untested, is that it was calibrated on a different ocean-heat driver.
  Tony would know.
- The fetch script documents all of this; the URL is still owed for OUR file.

## 2026-10-08d — Likelihood weights for the record-crossing draws: the record rejects every onset through 2025

**Asked (Marcus 10-08):** run the likelihood-weighting test on the crossing draws, in place of a hard cutoff at the
record end.

**Method** (`scripts/record_likelihood_weights.jl`):
- L27 was sampled with λ and T_crit at their paleo medians, and the propagated pair is attached afterwards. So the
  posterior conditioned on the record is reached exactly by importance weights,
  w = L(θ, λ_i, T_crit,i) / L(θ, λ_med, T_crit,med).
- Only the Antarctic level term can move: the SMB anchor is a 1979–2008 mean, before any onset. So
  log w = log L_ais(fast dynamics on) − log L_ais(T_crit = +∞).
- The likelihood is the calibrator's own (`hetero_logl_ar1`, L27's `--obs-corr-len=100`, target 1900–2025), with each
  draw's `sd_ais` and `rho_ais`.

**Gates (PASS):**
- [RESID]: the on-residuals equal the Table 5 residual file (gate-4 identical to the research code) for all 10,000
  draws.
- [W=1]: all 9,996 non-firing draws have identical on and off residuals, so w = 1 exactly.
- The comparator has power: the 4 firing draws do differ.

**Calibration climate (the posterior's own importance weights):**

| file row | onset | excess AIS 2025 | Δ log L | w |
|---|---|---|---|---|
| 9713 | 2021 | 1.33 cm | −802 | 0 |
| 4618 | 2025 | 0.63 cm | −397 | 6e-173 |
| 4588 | 2025 | 0.20 cm | −46.2 | 9e-21 |
| 4182 | 2026 | 0 | 0 | 1 exactly (after the record) |

- **Pre-check for 4588:** a 0.20 cm jump against an innovation sd of 0.023 cm is 8.6 sd, so Δ log L ≈ −37. The
  measured −46 agrees.
- **So the likelihood reproduces `record_end()` as the cutoff:** an onset by 2025 is rejected, and an onset in 2026 is
  untouched. `posterior(; drop_record_crossings=true)` at its default is the likelihood-weighted posterior, up to weights
  ≤ 1e-20.

**Joint arm (each pair on its own spliced config, the same ratio; a flagged choice):**
- All 3 SSP pairs and both vv pairs that fire by 2025 have Δ log L between −609 and −1703 (excess 1.1–2.3 cm).
- Every 2026 onset gets w = 1.

**Noise-model sensitivity (reported, not used):**
- The band collapses from 0.091 cm (2014–18) to the 0.05 cm floor at the 2019 GRACE-FO splice, and a correlated band
  (L = 100) cannot absorb a one-year jump. Three lenient models: diagonal band; diagonal band widened to its 2014–18
  level; diagonal band tripled.
- Under them, row 4588 gets w = 4e-4 / 0.05 / 0.35. Its rejection is therefore a property of L27's error model.
- Every other crossing, including all five joint pairs, is rejected under all three (w ≤ 2e-14).
- 4588 is not a projection row, so the paper's joint arm is the same under every variant.

## 2026-10-08c — v1.2-DEV: opt-in conditioning of the paleo draws on the record; v1.1.0 tagged

**Ruling (Marcus 10-08, decision 2 of SLR-RFF-BRICK handoff 2026-10-08):** no change to the paper's method or numbers.
Conditioning is offered to package users as an OPTION, in v1.2, not v1.1.

**Versioning:** `v1.1.0` tagged at 78b03c6, the commit certified against the research v1.1 outputs. `Project.toml` →
`1.2.0-DEV`. Defaults are unchanged, so every v1.1 product is unchanged. (`CITATION.cff` still says 1.0.0, a staleness
that predates this entry; left for the release.)

**Code:**
- `src/record_crossing.jl`:
  - `fastdyn_onset_year(row, gmst, years)`: the first year fast dynamics is on. It uses MimiBRICK's trigger,
    T_ant(t) = (GMST(t−1) − intercept) / coefficient > T_crit, with the arithmetic `apply_draw!` sets up.
  - `record_end()`: the Antarctic target's last observed year (2025), read from `calibration_targets.csv`.
  - `record_crossings(post; ...)`.
  - The record climate is `RECORD_FORCING = "ssp245harm"`, the calibration's own hindcast forcing.
- `posterior(; drop_record_crossings=false, record_last_year=nothing)`.
  - With the option on, the selected draws whose fast dynamics fires by the record end on that climate are deleted,
    and the dropped file rows are logged. That is rejection sampling, i.e. exact conditioning.
  - It drops file rows 4588, 4618 and 9713 (onsets 2025, 2025, 2021). With `record_last_year=2026` it also drops 4182.
    None is in `1:5:10000`.
- `scripts/record_crossings.jl` writes the onset year of all 2,000 projection draws, per scenario, in three climates:
  the record climate, the scenario mean (the fixed arm), and the draw's own spliced config (the joint arm).

**Finding: the joint arm DOES contain crossing draws.**
- The handoff's "none of the 4 is among the 2,000 projection rows, so no reported projection contains one" was
  measured on the mean climate only.
- On each draw's own FaIR config:
  - **SSPs:** 3 of the 2,000 fire by 2025 (in 2021–2023) and 7–8 by 2026.
  - **van Vuuren markers:** 2 fire by 2025 and 5 by 2026.
  - On the mean climate, none.
- Across all 10,000 draws paired the same way: 20 by 2025 for the SSPs, 8 for the markers.
- Effect on the reported numbers: SLR-RFF-BRICK CHANGELOG 2026-10-08e.

**Checks** (`test/conditioning/check_record_crossing.jl`, PASS):

| check | what it tests |
|---|---|
| [AGREE-RECORD] | `fastdyn_onset_year` equals the model's own trigger (first year `disintegration_rate ≠ 0`) on all 10,000 draws, 1850–2026: 0 mismatches. The crossing set is the documented one. |
| [AGREE-JOINT] | The same on the SSP2-4.5 joint arm, 2,000 draws on their own configs: 0 mismatches. |
| [BOUNDARY] | T_crit set to a draw's own record-maximum T_ant must not fire (the trigger is strict); 1 ulp below, it must fire, in the test and in the model, on 5 draws. |
| [INERT] | `drop_record_crossings=false` returns the v1.1 posterior; on `1:5:10000` the option drops nothing. |
| [DELETION] | The conditioned frame is the unconditioned one with exactly the crossing rows deleted, every kept row unchanged. |
| tag diff | Default `posterior()`, `rows=1:5:10000` and the v1.0 joint-arm call are `isequal` to the same calls at tag `v1.1.0`. |

**Mutation test (2/2 killed):**
- `>` → `>=` is caught only by [BOUNDARY] (5 failures). No real draw sits exactly on the threshold, so the agreement
  checks are blind to it. That is why [BOUNDARY] exists.
- The off-by-one year (`gmst[j]` for `gmst[j−1]`) is caught by every block (10 failures).
- Both mutants ran on scratch copies; `src/` was never edited for them.

## 2026-10-08b — v1.1 re-certified against the research v1.1 outputs: gates 2–4 PASS, plus a new panel gate

The research rerun is SLR-RFF-BRICK `run_ladrillo_v11_rerun_20261008.sh`, 0 failures. It passed three gates of its own:
- [V1-REGRESSION]: 16 shipped v1.0 files reproduced exactly;
- [CONTROL-EXACT]: on 15 SSP arms;
- Table 5 inputs: byte-identical.

New `test/identity/run_products.sh` produces every compared product with this package's defaults in three streams
(24 runs, about 9 min). New `compare_all.sh` runs the gates.

| gate | compared | result |
|---|---|---|
| 2 products | joint draws, Ladrillo + BRICK 2.0, 3 SSPs + 7 van Vuuren markers | PASS |
| 2b panel (new, `compare_panel.jl`) | `project_fixed.jl` vs the research fixed panel, 5,598 rows × 7 columns | PASS |
| 3 postpred | posterior predictive (moves under v1.1 only in 2021–2026, through the threshold-crossing draws) | PASS |
| 4 hindcast | BRICK 2.0 hindcast + IC residuals | PASS (the residuals did not move) |

- **Gate 5 (Table 5) is not re-run:** its inputs are byte-identical to v1.0 on the research side, and gate 4 ties this
  package's residuals to them.
- **Panel-gate power:** against the v1.0 panel it fails on 2,736 cells.
  - Antarctica and the total differ in every row.
  - Land water differs by at most 2.8e-14 cm, from re-referencing a shifted series.
  - Glaciers, Greenland and TE differ nowhere, which is the ordering again.
- **Tried and fixed:** the first `run_products.sh` used `$BASHPID`, which does not exist in macOS's bash 3.2. Every
  step failed instantly under `set -u`. Each stream now names its log explicitly.

## 2026-10-08 — v1.1: paleo single assignment and land-water frame fix (defaults); v1.0 tagged and kept reachable

**Rulings (Marcus 10-08):**
- The GMD paper moves to v1.1.
- The research drivers produce the numbers, and this package is re-certified against them.
- `v1.0.0` is tagged at a2bb3b0, certified against SLR-RFF-BRICK 5184a57.

**Code:**
- `Project.toml` 1.1.0.
- `PALEO_ASSIGNMENT = :single`: `posterior()` assigns over the whole file, then takes `rows`. `paleo = :v1` with
  `paleo_block` reproduces v1.0.
- `LWS_OBS_ANCHOR = :zero_at_first_year`, plus `LadrilloModel(...; lws_anchor=)`.
- The three projection scripts take `--paleo=`, `--lws-anchor=` and `--v1.0`. Non-default settings go into the
  output names and the provenance.

**Identity gate:**
- `cases.jl` now carries the paleo rule and the anchor per case. K1–K4 are the v1.0 regression; `_v11` runs the new
  defaults.
- `dump_old.jl` serves the cases that match the research kernel's settings, which it reads from the environment.

**Results:**
- **v1.0 regression PASS:** this package with its v1.0 options is bit-identical to the v1.0 research kernel on K1–K4,
  using freshly regenerated references (sha256 in `out/v10_research_5d3add0/SHA256SUMS`).
- **v1.1 PASS:** bit-identical to the edited research kernel on K1_v11–K4_v11. Two independent implementations of the
  single assignment agree.
- **`fix_effects.jl` PASS, mutation-tested (2/2 mutants killed).**
  - K1 → K1_v11 isolates the land-water fix.
    - Glaciers, Greenland and TE are identical.
    - Land water moves only before 1900, by exactly the 1.5646 cm step.
    - Antarctica moves at ~1e-4 of itself, within 3× of the analytic (1+γ)/depth × step: max 0.014 cm at 2100 and
      0.049 cm at 2300.
  - K3 → K3_v11 isolates the paleo fix on the hindcast.
    - One of the 1,000 draws moves (from 2024), and only through a fast-dynamics threshold crossing: with T_crit = +∞
      in both assignments it agrees exactly.
  - **Tried and abandoned:** the first version of this gate bounded Antarctica by the 10-01 "≤ 0.012 cm" and FAILED at
    0.049 cm.
    - That number was the transient observed-vs-stylized land-water swap. The step is a permanent offset, and a
      linearity test (doubling the step doubles ΔAIS, ratio 2.000) confirmed the mechanism.
    - Its hindcast block asserted "identical" and FAILED.
    - Both were replaced by the physics and channel gates above.

## 2026-10-07e — Table 5 ported; identity gate 5 PASS

`python/ic_ladrillo_vs_brick20.py` (Table 5) reproduces the research `ic_ladrillo_vs_brick20_L27{,_rho0.95,_rho0.9}`
outputs exactly: 3 × 182 table values and 3 × 420,000 per-draw likelihood values
(`test/identity/compare_ic_table.py`). It runs on 10,000 draws per model, about 40 min per rho bound with 4 workers.

**Comparator power:** a 1-ulp change to one table value, or to one per-draw value, is caught (exactly 1 difference
each).
- **Tried and fixed:** the first version parsed with pandas' default float parser. Re-serialising a file through that
  parser changed 310 other near-zero values (max 4e-25). The parser is not exact, so it could in principle merge two
  values 1 ulp apart. Both files are now read with `float_precision="round_trip"`. The Julia comparators use CSV.jl,
  which is exact.

`requirements.txt` pins numpy 2.4.2, pandas 2.3.3 and scipy 1.17.1, the environment the gate ran in. The README gains
the hindcast and Table 5 reproduction steps and the list of product gates.

## 2026-10-07d — BRICK 2.0 hindcast and IC residuals ported; shared targets and BRICK 2.0 code moved into src/

**New in `src/`:**
- `targets.jl`: `calibration_targets()`, the observations, the per-year sigma the likelihood uses, and the observed
  land water. One definition, used by every hindcast script.
- `brick20.jl`: `brick20_posterior`, `brick20_model`, `set_brick20_forcing!`, `apply_brick20!` and
  `brick20_sea_level`. `brick20_model` seeds the global RNG to BRICK20_SEED immediately before `get_model`.
  `scripts/brick20_joint.jl` now uses them, and its provenance column records the seed.

**New scripts:**
- `posterior_predictive_brick20.jl`: the BRICK 2.0 hindcast bands, on the same conventions as the Ladrillo check.
- `ic_hindcast_residuals.jl`: per-draw residuals of both models (10,000 draws each), the input to Table 5. Its built-in
  gate requires the residual medians to reproduce both posterior-predictive p50 series to 1e-9.
- `python/ic_ladrillo_vs_brick20.py`: Table 5 itself, ported (gate 5 pending).

**Regression after the refactor:**
- Gate 2: all 7 vv markers, both models; 504,000 Ladrillo and 252,000 BRICK 2.0 values identical.
- Gate 3: identical.

**Identity gate 4** (`test/identity/compare_hindcast.jl`): `postpred_brick20` and the two 10,000 × 644 residual files
are identical to the research `postpred_oldbrick`, `ic_hindcast_residuals_{ladrillo_L27,brick20}` and
`ic_hindcast_obs_sigma`.

**Mutation test of gate 4:**
- BRICK 2.0 arm re-referenced to 1995-2014: killed (6,270,000 values).
- Ladrillo land water `:observed` instead of `:central`: killed (2,510,000 values, up to 0.018 cm). Only AIS and the
  total move, through MimiBRICK's global-sea-level feedback onto Antarctica.
- **The planned no-power control was NOT one.** I expected the BRICK 2.0 get_model seed not to matter, because the
  hindcast never reads BRICK's own land water. It does matter: it moves 110,000 values by up to 8.8e-4 cm, in
  `ais_2020-2025` and `total_2020-2024` only. BRICK's random land-water realisation starts in 2019 and reaches
  Antarctica through the same sea-level feedback. So the seed is a real input of the BRICK 2.0 hindcast, and recording
  it in the artifact is required, not cosmetic.

**Also:** CITATION.cff (author Marcus C. Sarofim, NYU Marron Institute; the paper DOI is added when it exists).

## 2026-10-07c — posterior predictive ported; identity gate 3 PASS

**`scripts/posterior_predictive.jl`** ports the research `posterior_predictive_ladrillo.jl` on its L27 path. It runs the
posterior over 1850-2026 on the calibration forcing (ssp245harm), with land water `:central`, the tap off, and sea level
re-referenced to the calibration window 1995-2005. The parameter and predictive bands, bias and coverage go to
`outputs/hindcast/`.

**Kernel additions**, both off by default so projections are unchanged:
- `LadrilloModel(...; ref=)`: the window `sea_level` re-references to. The default is `REF_PERIOD`; the new constant
  `FIT_REF` = 1995-2005 is the calibration's window.
- `sea_level(...; funch=u)` and `funch_unit(years)`: the uncharted-ice hindcast term. It never enters the Mimi graph.

**Data:** `observations/calibration_targets_gsic_adj.csv` was imported (the glacier target with the RGI-19 seam removed).
The manifest's other rows are unchanged except the source commit (0c4de0a → c35513a; data identical).

**Identity gate 3** (`test/identity/compare_postpred.jl`) compares with the shipped research
`postpred_L27_{components_timeseries,bias,coverage}.csv` using `isequal` on every shared column. All 127 + 30 + 5 rows
are IDENTICAL. The package files add a `provenance` column to bias and coverage, which the research files lack.
- **Mutation test**, both mutants killed:
  - noise seed 2026 → 2027: 1,016 band values and 2 coverage values differ. Bias is unchanged, as it should be,
    because bias uses the parameter band only.
  - F_unch dropped from the total only: 381 + 12 + 2 values differ.
- **v1.0 simplifications:** the research driver's `--tag`, variant, AIS-ramp and `--ais-fit-from` recovery are dropped.
  L27 fits all four component series from 1900; which series are fitted is still read off the posterior's `sd_*`/`rho_*`
  columns.

## 2026-10-07 — package built from the research repository; certified bit-identical

**What this is.** A clean, L27-only `Ladrillo.jl` (Marcus 2026-10-07), built from `SLR-RFF-BRICK` (`ladrillo-dev`,
commit 0c4de0a). The glacier blocks are renamed SLOWP→SLOWG and FAST→FASTG; the rename is label-only
(`data/MANIFEST.csv`, `tools/import_l27_data.py`).

**Structural simplifications**, each covered by the identity gate below:
- Greenland is a dedicated two-basin component. It replaces the three-basin component run with a zero middle share.
- The shipped tap is now structural: a two-stage cascade with the whole-sheet clamp.
- The glacier regrowth divisor is dropped (it was 1.0).
- Every non-L27 variant path is removed.

**Identity gate 1 (kernel)** compares the package with the research kernel on every component × year × draw, using
`==`. All four cases pass: ssp585 on all 10,000 draws with the tap on; ssp245 with no tap; the hindcast under :central
land water; and the vvH joint-climate arm.

**Identity gate 2 (products)** compares the package's joint-climate drivers with the research drivers' outputs. On
vvVL, all 72,000 Ladrillo values and all 36,000 BRICK 2.0 values are identical. The research driver read the raw
2 GB chains; the package reads the 10k subsample. That proves the decision-2 switch moved nothing:
`subsample[1:5:end]` is bit-identical to the joint arm's chain draws.

**Mutation test** (case K4). Killed: a 1-ulp change to the active-basin share, a re-associated glacier sum, and a
paleo assignment in one block. Survived, as expected: two no-power controls.
- **Tried and replaced:** the first reassociation mutant reordered the Greenland headroom sum. It survived because the
  whole-sheet clamp never binds, so it is kept as a no-power control.

**The MimiBRICK build scenario was checked:** building with ssp585 vs ssp245 changes 0 of 200 draws × 6 components ×
451 years. So every product is built on ssp245, which also lets the van Vuuren scenarios run.

**Decisions (Marcus 2026-10-07):**
- **Paleo attachment.** Keep v1.0's rules: the joint arm assigns per chain (`paleo_block = 500`); the SSP panel
  assigns once over its 2,000 rows. Move to a single assignment at v1.1.
  - Measured on ssp585: medians move ≤ 1.3 bootstrap se.
  - The 2100 p95 moves −3.6 cm (total) and −5.0 cm (AIS), about 2.5 se.
  - Tool: `scripts/project_joint.jl --paleo=single`.
- **BRICK 2.0 posterior.** Not redistributed: `tools/fetch_brick20_posterior.sh` fetches it and checks its sha256,
  but the URL is still to come. The file was purged from this repository's local history before any remote existed.
- **The van Vuuren forcing is the CMIP7 basis** (the research repo's 2026-10-02 regeneration). The research repo's vv
  arms and FACTS are being re-run on it.
