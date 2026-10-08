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

**This is version 1.0, the posterior of the model-description paper (calibration L27).** The posterior is
`data/posterior/ladrillo_L27_posterior.csv`: 10,000 draws pooled from four chains.

## Install

Julia 1.12 (tested on 1.12.6).

```bash
git clone <this repository> Ladrillo && cd Ladrillo
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

*In progress (v1.0 release): the calibration, the FaIR forcing builders and the figure scripts.*

## Validation

`test/identity/` holds the gate that certified this package against the research code that produced the paper:
every component, every year and every draw must be bit-identical (`==`), on the full 10,000-draw posterior and on
the joint-climate arm. The gate is mutation-tested (`test/identity/mutation_test.sh`, result in
`test/identity/mutation_result_K4.txt`): a one-ulp change to a Greenland volume share, a re-associated glacier sum,
and a different paleo-draw assignment are each caught; two mutants on paths that never bind in v1.0 (the glacier
floor, the Greenland whole-sheet clamp) survive, as they should.

The paper's products are gated the same way, against the research drivers' own output files:
- `compare_products.jl`: the joint-climate projections on all seven van Vuuren scenarios, both models.
- `compare_postpred.jl`: the posterior predictive.
- `compare_hindcast.jl`: the BRICK 2.0 hindcast and the residuals behind Table 5.
- `compare_ic_table.py`: Table 5 itself.

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
python/              the Table 5 information-criterion test
test/                tests and the identity gate
tools/               one-time import of the data from the research repository
```

## Known issues

- **Land-water frame step (fix scheduled for v1.1).** The observed land-water series is in the calibration target's
  1995-2005 frame, while the other components start at 0 in 1850, so the land-water component steps by +1.565 cm in
  1900 and Antarctica's sea-level feedback sees the step. The effect is at most 0.012 cm of Antarctic contribution at
  any year to 2300, and reported totals are re-referenced. v1.0 keeps the step so the paper's numbers reproduce
  (`Ladrillo.LWS_OBS_ANCHOR`).

## License

Code: MIT (`LICENSE`).
