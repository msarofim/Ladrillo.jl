# Changelog

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
