# Changelog

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
