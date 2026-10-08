# Ladrillo.jl

Ladrillo is a calibrated, reduced-complexity emulator of global mean sea-level rise. It is
[MimiBRICK](https://github.com/raddleverse/MimiBRICK.jl) v2.0.0 with two components replaced and the whole model
recalibrated jointly against observed component budgets:

| component | in Ladrillo |
|---|---|
| glaciers and ice caps | three reservoirs (R19 Antarctic periphery; SLOWG, the slow-responding RGI regions 03/09/07/06; FASTG, the other 13 regions), each on its own area-weighted temperature driver |
| Greenland | two basins (active SW/CW/CE/SE/NW and high NO/NE, Mouginot et al. 2019 sectors), fast (surface mass balance) and slow (dynamic) channels, plus a prior-specified above-threshold discharge channel |
| Antarctica | MimiBRICK's DAIS component; calibrated with paleo priors; fast dynamics propagated from the paleo ensemble |
| thermal expansion, land water | MimiBRICK's; land water follows the observed record (Frederikse et al. 2020; GRACE/GRACE-FO) |

The model is driven by prescribed global mean surface temperature and ocean heat content. The repository ships
FaIR 2.2.4 (calib 1.6.0) forcing for SSP1-2.6, SSP2-4.5, SSP5-8.5 and the seven van Vuuren et al. (2025) scenario
markers, as both the 841-config ensemble mean and the full ensemble.

**This is version 1.1, on the posterior of the model-description paper (calibration L27).** The posterior is
`data/posterior/ladrillo_L27_posterior.csv`: 10,000 draws pooled from four chains. Version 1.1 changes no parameter; it
fixes two things in how projections use the posterior (see *Versions* below). Version 1.0, the code as submitted with
the paper, is the git tag `v1.0.0`, and both of its behaviours remain reachable as options.

## Install

Julia 1.12 (tested on 1.12.6).

```bash
git clone https://github.com/msarofim/Ladrillo.jl.git Ladrillo && cd Ladrillo
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` pins every dependency, MimiBRICK to the `v2.0.0` tag.

## Quick start

```julia
using Ladrillo
gmst, ohc = mean_forcing("ssp245")          # FaIR ensemble mean, 1850-2300
lm = LadrilloModel(gmst, ohc)               # land water :observed, Greenland discharge channel on
post = posterior()                          # 10,000 draws
run_draw!(lm, post[1, :])
sea_level(lm, :total)                       # cm rel. 1995-2014, one value per year 1850-2300
```

Components: `:glaciers`, `:gis`, `:ais`, `:te`, `:lws`, `:total`, and the glacier reservoirs `:gsic_r19`,
`:gsic_slowg`, `:gsic_fastg`. Any forcing on any year span can be passed in place of `mean_forcing`.

## Reproducing the paper's projections

| product | command | runtime (one core, Apple M4) |
|---|---|---|
| projection with FaIR climate uncertainty (one scenario) | `julia --project=. scripts/project_joint.jl --scenario=ssp245` | ~5 min |
| BRICK 2.0 on the same climate | `julia --project=. scripts/brick20_joint.jl --scenario=ssp245` | ~1 min |
| the fixed-climate SSP panel | `julia --project=. scripts/project_fixed.jl` | |

Outputs go to `outputs/projections/`; every file carries a `provenance` column. Scenarios: `ssp126 ssp245 ssp585
vvVL vvLN vvL vvML vvM vvHL vvH`.

## Reproducing the hindcast comparison (Table 5)

Run in this order; each later step checks itself against the earlier ones.

| product | command | runtime |
|---|---|---|
| Ladrillo posterior predictive (bands, bias, coverage) | `julia --project=. scripts/posterior_predictive.jl` | ~2 min |
| BRICK 2.0 hindcast on the same conventions | `julia --project=. scripts/posterior_predictive_brick20.jl` | ~2 min |
| per-draw residuals of both models, 10,000 draws each | `julia --project=. scripts/ic_hindcast_residuals.jl` | ~5 min |
| Table 5, the information-criterion test (rho bound 0.99; 0.95 and 0.9 for the sensitivity) | `python python/ic_ladrillo_vs_brick20.py --rho-max=0.99` | |

The Python step needs `pip install -r requirements.txt` (NumPy, pandas, SciPy). `IC_WORKERS` sets its process count.

Outputs go to `outputs/hindcast/` and `outputs/tables/`. Runtimes were measured on an Apple M4 that was busy with other
jobs; expect them to be shorter.

Tables A1 and A2 (seconds, NumPy and pandas, to `outputs/tables/`): `python python/tables/ladrillo_prior_posterior_table.py` for A1; `python python/tables/diag_ais_block_pca.py` then `python python/tables/ladrillo_table_a2.py` for A2.

## Rebuilding the FaIR forcing

`tools/forcing/` rebuilds every file in `data/forcing/` from its inputs, in about a minute, and checks them
byte-for-byte against the shipped files (`check_forcing_identity.py`). `tools/forcing/README.md` covers the recipe,
the inputs, their checksums and their sources.
- ⚠ Two inputs are not public: the IIASA ScenarioMIP prerelease emissions, needed for the van Vuuren forcing, and an
  AR6 SSP2-4.5 file, which shapes the calibration forcing after 2100 only.

*In progress: the calibration code and the figure scripts.*

## Validation

**Run the tests** (from a fresh clone; nothing outside this repository is needed):

```bash
julia --project=. test/runtests.jl
```

That takes about 2 minutes on an Apple M4, and `julia --project=. -e 'using Pkg; Pkg.test()'` runs the same tests in
Pkg's sandbox. The suite re-runs 33 posterior draws in both projection arms (SSP2-4.5, SSP5-8.5), the full posterior
predictive and a BRICK 2.0 sample. It compares them with values the RESEARCH code wrote (`test/reference/`, sources and
sha256 in `test/reference/SOURCES.csv`). It also checks every data file against `data/MANIFEST.csv`, the component-sum
and splice identities, the v1.0 options, the record-crossing test, and Tables A1 / A2 (when Python with
NumPy and pandas is available).
- **On the reference runtime** (macOS arm64, Julia 1.12.6, default bounds checking) every comparison is exact (`==`).
  A one-ulp change to a live parameter fails it.
- **Elsewhere, including `Pkg.test()`, which forces `--check-bounds=yes`,** values must agree to a relative 1e-12. The
  last bits can move with code generation: 1.7e-15 measured under bounds checking.
- **BRICK 2.0** is compared only once its posterior is in place (see *BRICK 2.0 posterior* below), and so is that
  file's checksum. Until then a fresh clone reports those two checks as skipped (Julia's Test prints them as "Broken");
  a fetched file whose checksum is wrong fails.

`test/identity/` holds the gate that certified this package against the research code that produced the paper:
every component, every year and every draw must be bit-identical (`==`), on the full 10,000-draw posterior and on
the joint-climate arm, for both the v1.0 settings and the v1.1 defaults (`cases.jl`). The gate is mutation-tested (`test/identity/mutation_test.sh`, result in
`test/identity/mutation_result_K4.txt`): a one-ulp change to a Greenland volume share, a re-associated glacier sum,
and a different paleo-draw assignment are each caught; two mutants on paths that never bind in v1.0 (the glacier
floor, the Greenland whole-sheet clamp) survive, as they should.

The paper's products are gated the same way, against the research drivers' own output files:
- `compare_products.jl`: the joint-climate projections on all seven van Vuuren scenarios, both models.
- `compare_postpred.jl`: the posterior predictive.
- `compare_hindcast.jl`: the BRICK 2.0 hindcast and the residuals behind Table 5.
- `compare_ic_table.py`: Table 5 itself.

`test/conditioning/check_record_crossing.jl` checks the v1.2 record-crossing test (*Versions*).

## BRICK 2.0 posterior

The BRICK 2.0 comparison uses Tony Wong's post-PR#93 MimiBRICK posterior subsample (10,000 draws, sha256 `78247db8…`).
It is not redistributed here, and `tools/fetch_brick20_posterior.sh` still needs its public URL.
- ⚠ The public MimiBRICK v2.0.0 Zenodo record (10.5281/zenodo.20592337) ships a different calibration under the same
  file name. On this package's forcing its thermal-expansion hindcast is 4× worse.
- It is therefore not a substitute: the paper's BRICK 2.0 numbers rest on the sha256 above.

## Repository layout

```
src/                 the model: components/, land water, posterior, forcing helpers
data/posterior/      the v1.0 posterior
data/drivers/        observed glacier-region and Greenland temperature drivers (HadCRUT5-based)
data/constants/      fixed constants, the Greenland amplification shape, paleo draws
data/observations/   the calibration targets (and the glacier target with the RGI-19 seam removed)
data/forcing/        FaIR 2.2.4 (calib 1.6.0) GMST and ocean heat content
data/brick20/        the BRICK 2.0 comparison posterior
data/MANIFEST.csv    source file, commit and sha256 of every data file
scripts/             Julia drivers that produce the paper's products
python/              the Table 5 information-criterion test; tables/ (Tables A1, A2)
test/                runtests.jl + reference/ (self-contained), identity/ (the certification against the research repo)
tools/               data import, test-reference builder, BRICK 2.0 fetch; forcing/ (the FaIR builders)
```

## Versions

**v1.1** (2026-10-08). Same posterior; two projection-side fixes.
- **Paleo assignment** (`Ladrillo.PALEO_ASSIGNMENT = :single`). The Antarctic fast-dynamics pair (lambda, T_crit) is
  attached to each posterior draw as a joint draw from the paleo ensemble.
  - v1.0 made that assignment after selecting draws, so its joint-climate arm, which assigned chain by chain, gave all
    four chains the same 500 paleo rows.
  - v1.1 assigns once over all 10,000 draws, and any subset keeps the rows its draws have in the full table.
  - Measured on SSP5-8.5 (joint arm, 2026-10-07): medians move by at most 1.3 bootstrap standard errors; the 2100 95th percentile moves
    -3.6 cm (total) and -5.0 cm (Antarctica).
- **Land-water frame step** (`Ladrillo.LWS_OBS_ANCHOR = :zero_at_first_year`).
  - The observed land-water series is in the calibration target's 1995-2005 frame, while the other components start
    at 0 in 1850. v1.0 started the land-water component at the series' 1900 value, +1.565 cm, and Antarctica's
    sea-level feedback saw that step. v1.1 starts it at 0.
  - Glaciers, Greenland and thermal expansion are bit-identical.
  - Antarctica moves by about -1e-4 of its own contribution: at most 0.014 cm at 2100 and 0.049 cm at 2300 under
    SSP5-8.5.
- `--v1.0` on `scripts/project_joint.jl`, `project_fixed.jl` and `brick20_joint.jl` (or `posterior(paleo=:v1, ...)`
  and `LadrilloModel(...; lws_anchor=:v1_step)`) reproduces v1.0 bit for bit. The identity gate checks this.
- `test/identity/fix_effects.jl` separates the two fixes and checks each against what it should and should not move.

**v1.2** (in development, `1.2.0-DEV`). Adds an OPT-IN conditioning of the propagated paleo draws on the record; the
defaults, and every v1.1 product, are unchanged. `v1.1.0` is tagged.
- `posterior(; drop_record_crossings=true)` drops the draws whose Antarctic fast dynamics fires by the end of the
  Antarctic record (`record_end()`, 2025) on the calibration's own climate. That is rejection sampling, i.e. exact
  conditioning of the joint (posterior x paleo) draws on "no fast dynamics through 2025". It drops 3 of the 10,000 draws
  (`record_last_year=2026`: 4), none of them among the 2,000 projection rows.
- That tests ONE climate. In the joint-climate arm each draw runs on its own FaIR config, and more draws cross there
  (see *Known issues*). Condition such a run on each pair's own GMST with `fastdyn_onset_year(row, gmst, years)`, after
  pairing: the arm pairs draws with configs by position.
- `fastdyn_onset_year` is the model's own trigger with the model's arithmetic. `test/conditioning/check_record_crossing.jl`
  checks it against model runs (10,000 draws on the record climate, 2,000 joint draws) and at the threshold to 1 ulp,
  and is mutation-tested. `scripts/record_crossings.jl` lists the onset year of every projection draw in each climate.

## Known issues

- **The fast-dynamics threshold is crossed inside the hindcast by a few propagated draws.**
  - The calibration holds lambda and T_crit at their paleo medians, which never cross, but the paleo ensemble's
    low-T_crit tail does: 4 of the 10,000 draws (T_crit -16.4 to -16.9 C, DAIS scale) start fast dynamics in
    2021-2026 on the calibration climate, adding up to 1.6 cm of Antarctic sea level by 2026.
  - On its own FaIR config a draw can cross where it does not on the mean. In the joint-climate arm, 3 of the 2,000
    projection draws fire by 2025 under each SSP (2 under each van Vuuren marker), and 7-8 (5) by 2026; on the mean
    climate, none.
  - Hindcast products that subsample the posterior can therefore depend slightly on the paleo assignment.
  - v1.1 (and the paper) does not condition on the record; v1.2 offers it as an option (*Versions*).

## License

Code: MIT (`LICENSE`).
