# Changelog

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
