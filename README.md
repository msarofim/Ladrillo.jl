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

**This is version 1.0, the release the model-description paper reports, on its posterior (calibration L27).**
The posterior is `data/posterior/ladrillo_L27_posterior.csv`: 10,000 draws pooled from four chains. No version changes a
parameter.
- Versions 0.1 to 0.4 were internal builds, git tags `v0.1.0` to `v0.4.0` (tagged `v1.0.0` to `v1.3.0` before the
  release). See *Versions* below.
- Version 1.0 has v0.4's model and numbers. It adds the scripts for the paper's figures and renames the options that
  reproduce v0.1 (`--v0.1`, `paleo=:per_selection`, `lws_anchor=:obs_at_first_year`). The behaviours of v0.1 and v0.2
  remain reachable as options.

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

## Reproducing the hindcast comparison (Tables 4 and 5)

Run in this order; each later step checks itself against the earlier ones.

| product | command | runtime |
|---|---|---|
| Ladrillo posterior predictive (bands, bias, coverage) | `julia --project=. scripts/posterior_predictive.jl` | ~2 min |
| BRICK 2.0 hindcast on the same conventions | `julia --project=. scripts/posterior_predictive_brick20.jl` | ~2 min |
| Table 4, the hindcast scorecard (bias, RMSE, coverage by component and window) | `python python/tables/ladrillo_vs_brick20_scorecard.py` | seconds |
| per-draw residuals of both models, 10,000 draws each | `julia --project=. scripts/ic_hindcast_residuals.jl` | ~5 min |
| Table 5, the information-criterion test (rho bound 0.99; 0.95 and 0.9 for the sensitivity) | `python python/ic_ladrillo_vs_brick20.py --rho-max=0.99` | |

The Python step needs `pip install -r requirements.txt` (NumPy, pandas, SciPy). `IC_WORKERS` sets its process count.

Outputs go to `outputs/hindcast/` and `outputs/tables/`. Runtimes were measured on an Apple M4 that was busy with other
jobs; expect them to be shorter.

Tables A1 and A2 (seconds, NumPy and pandas, to `outputs/tables/`): `python python/tables/ladrillo_prior_posterior_table.py` for A1; `python python/tables/diag_ais_block_pca.py` then `python python/tables/ladrillo_table_a2.py` for A2.

`python python/tables/check_tables_identity.py <SLR-RFF-BRICK checkout>` compares Tables 4, A1 and A2 with the
research outputs behind the paper (every number exact on the versions in `requirements.txt`). Table A1's markdown shows
the four parameters sampled as log₁₀ in natural logs (footnote ᶜ); its CSV keeps the log₁₀ values.

## Rebuilding the FaIR forcing

`tools/forcing/` rebuilds every file in `data/forcing/` from its inputs, in about a minute, and checks them
byte-for-byte against the shipped files (`check_forcing_identity.py`). `tools/forcing/README.md` covers the recipe,
the inputs, their checksums and their sources.
- The van Vuuren forcing is built from the published Zenodo 20713982 v1.1.1 emissions.
- ⚠ One input is not public: an AR6 SSP2-4.5 file, which shapes the calibration forcing after 2100 only. Without it
  the build makes the other 41 of the 43 files and says which two it skipped.

## Rebuilding the calibration targets and drivers from their sources

`tools/targets/` rebuilds the eight observation-based files in `data/` (the calibration targets and their
seam-adjusted glacier version, the glacier and Greenland temperature drivers, the glacier reservoir constants, the
Greenland amplification prior and its CMIP6 shape) from the source datasets, in about 6 seconds, and checks them byte for
byte against the shipped files. `tools/targets/README.md` has the chain, the gate, its mutation test, what the port
changed, and the full sources table.

```bash
pip install -r tools/targets/requirements-targets.txt
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs             # what to fetch, from where
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --download  # the files with stable public URLs
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --extract   # needs the raw Berkeley Earth, GISTEMP, GTN-G files
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --verify
python tools/targets/build_targets.py --inputs-dir ~/ladrillo_inputs --work-dir /tmp/targets_build [--from-raw]
python tools/targets/check_targets_identity.py --work-dir /tmp/targets_build
```

- Small sources, and small extractions of large gridded products, are shipped in `data/sources/` (2.4 MB). Sources whose
  licence leaves redistribution in doubt are not shipped: `fetch_inputs.py` says where to get them, and
  `tools/targets/inputs.sha256` pins each file, so the builder refuses any other version.
- With `--from-raw` the chain starts from the raw gridded products (HadCRUT5, Berkeley Earth, GISTEMP, the GRACE mascon
  grid, the Dangendorf fields) and re-makes every shipped extraction byte for byte first.
- Several sources are rolling products that their publishers overwrite (NOAA STAR, NOAA NCEI, Berkeley Earth, GISTEMP).
  The exact files these targets were built from are not served any more; the pins name them.
- No builder exists for four constants (the central row and the three paleo files); they are frozen data with provenance.

<!-- SOURCES TABLE (compact): generated by `python tools/targets/target_inputs.py --write-readme`; do not edit by hand -->
| source | kind | licence | feeds |
|---|---|---|---|
| Frederikse et al. (2020): sea-level budget: global basin time series (xlsx) and the 5000-member weighted ensemble (nc) | fetch yourself | original record 3862995: 'Other (Not Open)', all rights reserved; record 10621070: CC BY 4.0, which cannot relicense Frederikse's files | every component target (ais, gis, gsic, steric, lws) 1900-2018 and its sigma; the total's budget-closure sigma; the land-water splice offset |
| Dangendorf et al. (2024), corrected global file: Kalman-smoother GMSL reconstruction, global series with standard errors (v2) | shipped | none stated; shipped by Marcus's ruling (2026-11-12) | dang_sig, dang_lo, dang_hi 1900-2021 (the total target's own SE) |
| Dangendorf et al. (2024): annual GMSL = cos(lat)-weighted mean of the HR field (extraction) | shipped | CC BY 4.0 (Zenodo API) | dang (the total target) 1900-2021 |
| JPL GRACE/GRACE-FO mascons RL06.3Mv04: Antarctic and Greenland mass time series | shipped | NASA EOSDIS data-use guidance: openly shared without restriction (CC0 unless marked) | ais and gis 2019-2026; the mask gates of the land-water extraction |
| JPL GRACE/GRACE-FO mascons RL06.3Mv04: monthly mass over global land minus both ice sheets (extraction) | shipped | as above (NASA, without restriction) | lws 2019-2023, held flat to 2026 |
| GlaMBIE Team (2025): regional calendar-year glacier mass change (results, calendar_years/ of glambie_data.zip) | shipped | WGMS data policy: open access with citation (= CC BY 4.0) | gsic 2019-2023 (global minus region 5); glacier block weights, modern rates, 2020 losses; the region-19 seam; glaciers removed from the land-water extension |
| NOAA NCEI: World Ocean 0-2000 m thermosteric sea level, yearly | fetch yourself (rolling) | none stated on the pages read | steric 2019-2025 |
| NOAA STAR Laboratory for Satellite Altimetry: global mean sea level, seasonal signals removed; annual means of the 2026-05 release | fetch yourself (rolling) | no licence; LSA data-use policy (2008) requires a credit line | dang 2022-2024 |
| GlacierMIP3 (Zekollari et al. 2025): Tables S1a and S3, region characteristics | shipped | CC BY 4.0 (Zenodo API) | glacier inventory, response times, ISIMIP3 amplification |
| GlacierMIP3 (Zekollari et al. 2025): committed-loss ladder per Ladrillo glacier reservoir (extraction) | shipped | CC BY 4.0 (Zenodo API) | the glacier equilibrium curves (b, T_off) and their committed rungs |
| Met Office HadCRUT5.0.2.0 (Morice et al. 2021): glacier-area-weighted annual temperature per RGI region, and the global mean (extraction) | shipped | Open Government Licence v3 (Crown copyright); acknowledgement line required | the three glacier reservoir drivers; the observational amplification fit |
| Met Office HadCRUT5.0.2.0: Greenland land-zone and global annual temperature (extraction) | shipped | Open Government Licence v3 (the zone weights use the Berkeley Earth 1-degree land fraction and the GTN-G region-05 polygon; no Berkeley values are in the file) | the Greenland driver; the amplification prior and anchor |
| Berkeley Earth: Greenland land-zone and global annual temperature (extraction, made by fetch_inputs.py --extract) | fetch yourself (rolling) | CC BY-NC 4.0 (berkeleyearth.org/data) | the amplification prior and anchor (one of three products) |
| NASA GISS GISTEMP v4: Greenland land-zone and global annual temperature (extraction, made by fetch_inputs.py --extract) | fetch yourself (rolling) | none stated (citation requested) | the amplification prior and anchor (one of three products) |
| CMIP6 (40 models): Amon tas, annual global and Greenland-zone means, historical + ssp126/245/585 (extraction) | shipped | CC BY 4.0 for all groups since 2022-10 (CMIP6 Terms of Use v6.2); GISS sources CC0 | the Greenland amplification shape S(dT) |
<!-- END SOURCES TABLE (compact) -->

## Reproducing the calibration

`scripts/calibrate.jl` samples the L27 calibration, one chain per call:

```bash
julia --project=. --threads=1 scripts/calibrate.jl 2000000 2026
```

The production run is four chains (seeds 2026-2029) of 2,000,000 iterations. Each takes about 3 h on an Apple M4 and
writes a ~2 GB chain to `outputs/calibration/`, with the adapted proposal covariance and a provenance file.
- `src/calibration.jl` holds the log-posterior, the parameter table, the start rows and the proposal seed. It ports the
  research calibrator (`julia/calibrate_mcmc_ext.jl` with the flags of `run_mcmc_L27.sh`) onto this package's model.
- `data/calibration/` holds its three frozen inputs: the four start rows (posterior draws of a tuning run), the proposal
  covariance adapted by the L26 run, and the Greenland amplification prior.
- `test/calibration/gate_chain_prefix.jl` (about 2 min) checks it against the production chains: the log-posterior at
  the four start rows, and the first 20,000 rows of the seed-2026 chain (4,603 accepted moves), are bit-identical to
  production. `test/calibration/mutation_test.sh` checks that the gate catches a changed sampler target, a tighter
  bound on the AR(1) rho and two changed likelihood constants.
- Bit-identity holds on the reference runtime only (Apple M4, Julia 1.12.6, the Manifest's versions, one BLAS thread).
  Elsewhere a single changed accept decision separates the chains, and the agreement is statistical.

`scripts/postprocess_calibration.jl` turns the four chains into the posterior:

```bash
julia --project=. --threads=1 scripts/postprocess_calibration.jl --accept-slr
```

It keeps the second half of each chain and checks that no chain is in the Greenland noise mode. It then computes R-hat
and ESS for every parameter, and for projected SSP2-4.5 sea level at 2100 and 2150 (400 draws per chain). The posterior
is every 400th pooled draw (10,000 draws, no random numbers), with the format of `data/posterior/`. L27 was accepted on
projected sea level (R-hat 1.001 and 1.002). Eight Antarctic parameters did not pass the parameter test (R-hat up to
1.31, ESS 20-244): five of the seven geometry parameters, antarctic_alpha, anto_alpha and anto_beta. `--accept-slr`
accepts on the projection in that case.
Without the flag the script writes nothing unless every parameter passes; `--force` writes with a `_NOTCONVERGED`
suffix. The production run takes about 25 minutes on an Apple M4 (21 of them the parameter ESS), 3.3 GB peak memory.
- `test/calibration/gate_postprocess.jl <SLR-RFF-BRICK repo>` runs it on the production chains (sha256 checked) and
  requires the posterior to equal `data/posterior/ladrillo_L27_posterior.csv` byte for byte. The sea-level R-hat table,
  the pooled covariance and the logged diagnostics (noise check, per-chain projected quantiles, the eight parameters'
  R-hat and ESS) must equal the research run of 2026-09-20. `test/calibration/mutation_test_postprocess.sh` checks that
  each of the four comparisons catches a one-line change to its own step.
- The calibration targets and drivers it reads are rebuilt by `tools/targets/` (*Rebuilding the calibration targets
  and drivers from their sources*).

## Reproducing the figures (Figs. 1–6)

Run these after the projection products for the seven van Vuuren scenarios (`project_joint.jl` and `brick20_joint.jl`,
above) and the two posterior-predictive steps (Fig. 1). Each takes a few seconds and needs matplotlib
(`pip install -r requirements.txt`).

| figure | command |
|---|---|
| Fig. 1, hindcast by component | `python python/figures/plot_hindcast_components.py` |
| the four-model comparison table that Figs. 2, 3 and 6 read | `python python/figures/vv_model_comparison.py` |
| Figs. 2 and 3, the four models at 2100 and 2300 (and 2150) | `python python/figures/plot_model_comparison_components.py` |
| Fig. 4, components to 2300 | `python python/figures/plot_future_components.py` |
| Fig. 5, glaciers | `python python/figures/plot_vv_gsic_ladrillo.py` |
| Fig. 6, response to the High minus Very Low scenario | `python python/figures/plot_vv_responsiveness.py` |

Figures go to `outputs/figures/`, each with its caption text in a `.caption.txt` beside it, and tables to
`outputs/tables/`. Figs. 2, 3 and 6 also show FACTS and MAGICC-SLR, which this repository does not run; their outputs
ship as data in `data/comparison/` (see its README).

`python python/figures/check_figures_identity.py <SLR-RFF-BRICK checkout>` compares every PNG byte for byte, every
caption and both tables with the research renders in the paper. With the versions in `requirements.txt` all 16 checks
pass; other matplotlib versions can move pixels without any change in the data.

## Validation

**Run the tests** (from a fresh clone; nothing outside this repository is needed):

```bash
julia --project=. test/runtests.jl
```

That takes about 2 minutes on an Apple M4, and `julia --project=. -e 'using Pkg; Pkg.test()'` runs the same tests in
Pkg's sandbox. The suite re-runs 33 posterior draws in both projection arms (SSP2-4.5, SSP5-8.5), the full posterior
predictive and a BRICK 2.0 sample. It compares them with values the RESEARCH code wrote (`test/reference/`, sources and
sha256 in `test/reference/SOURCES.csv`). It also checks every data file against `data/MANIFEST.csv`, the component-sum
and splice identities, the v0.1 options, the record-crossing test, the calibration's log-posterior at the four chain
starts and its parameter table, and Tables A1 / A2 (when Python with NumPy and pandas is available).
- **On the reference runtime** (macOS arm64, Julia 1.12.6, default bounds checking) every comparison is exact (`==`).
  A one-ulp change to a live parameter fails it.
- **Elsewhere, including `Pkg.test()`, which forces `--check-bounds=yes`,** values must agree to a relative 1e-12. The
  last bits can move with code generation: 1.7e-15 measured under bounds checking.
- **BRICK 2.0** is compared only once its posterior is in place (see *BRICK 2.0 posterior* below), and so is that
  file's checksum. Until then a fresh clone reports those two checks as skipped (Julia's Test prints them as "Broken");
  a fetched file whose checksum is wrong fails.

`test/identity/` holds the gate that certified this package against the research code that produced the paper:
every component, every year and every draw must be bit-identical (`==`), on the full 10,000-draw posterior and on
the joint-climate arm, for both the v0.1 settings and the v0.2 defaults (`cases.jl`). The gate is mutation-tested (`test/identity/mutation_test.sh`, result in
`test/identity/mutation_result_K4.txt`): a one-ulp change to a Greenland volume share, a re-associated glacier sum,
and a different paleo-draw assignment are each caught; two mutants on paths that never bind in v0.1 (the glacier
floor, the Greenland whole-sheet clamp) survive, as they should.

The paper's products are gated the same way, against the research drivers' own output files:
- `compare_products.jl`: the joint-climate projections on all seven van Vuuren scenarios, both models.
- `compare_postpred.jl`: the posterior predictive.
- `compare_hindcast.jl`: the BRICK 2.0 hindcast and the residuals behind Table 5.
- `compare_ic_table.py`: Table 5 itself.

`test/conditioning/check_record_crossing.jl` checks the v0.3 record-crossing test (*Versions*).

## BRICK 2.0 posterior

The BRICK 2.0 comparison uses the BRICK parameters of the joint SNEASY-BRICK posterior subsample published with
MimiBRICK v2.0.0 (Zenodo 10.5281/zenodo.20592337, `parameters_subsample_sneasybrick.csv`, 10,000 draws, CC-BY-4.0).
`tools/fetch_brick20_posterior.sh` downloads it to `data/brick20/` and checks its sha256 (`d74577f6…`); the tests then
run the BRICK 2.0 comparison.
- The same record also ships `parameters_subsample_brick.csv`, the standalone BRICK calibration. It is not this arm's
  file: on this package's FaIR forcing its thermal expansion is about 3× too small.
- Through v0.3.0 the arm used an earlier joint subsample from Tony Wong (sha256 `78247db8…`). The published file gives
  an SSP2-4.5 BRICK 2.0 total about 4 % lower at 2100 and 2300.

## Repository layout

```
src/                 the model: components/, land water, posterior, forcing helpers
data/posterior/      the L27 posterior (unchanged since v0.1)
data/drivers/        observed glacier-region and Greenland temperature drivers (HadCRUT5-based)
data/constants/      fixed constants, the Greenland amplification shape, paleo draws
data/observations/   the calibration targets (and the glacier target with the RGI-19 seam removed)
data/forcing/        FaIR 2.2.4 (calib 1.6.0) GMST and ocean heat content
data/brick20/        the BRICK 2.0 comparison posterior: fetched from Zenodo, not shipped (see *BRICK 2.0 posterior*)
data/comparison/     FACTS and MAGICC-SLR outputs for Figs. 2, 3 and 6 (static; see its README)
data/MANIFEST.csv    source file, commit and sha256 of every data file
data/calibration/    the calibration's start rows, proposal-covariance seed and Greenland amplification prior
data/sources/        the shipped inputs of tools/targets/ (sources and extractions; sha256 in tools/targets/inputs.sha256)
scripts/             Julia drivers that produce the paper's products, the calibration sampler and post-processing
python/              the Table 5 information-criterion test; tables/ (Tables 4, A1, A2); figures/ (Figs. 1-6)
test/                runtests.jl + reference/ (self-contained), identity/ (the certification against the research repo),
                     calibration/ (the sampler and post-processing gates)
tools/               data import, test-reference builder, BRICK 2.0 fetch; forcing/ (the FaIR builders);
                     targets/ (the calibration-target and driver builders, their inputs and gate)
```

## Versions

**v0.2** (2026-10-08). Same posterior; two projection-side fixes.
- **Paleo assignment** (`Ladrillo.PALEO_ASSIGNMENT = :single`). The Antarctic fast-dynamics pair (lambda, T_crit) is
  attached to each posterior draw as a joint draw from the paleo ensemble.
  - v0.1 made that assignment after selecting draws, so its joint-climate arm, which assigned chain by chain, gave all
    four chains the same 500 paleo rows.
  - v0.2 assigns once over all 10,000 draws, and any subset keeps the rows its draws have in the full table.
  - Measured on SSP5-8.5 (joint arm, 2026-10-07): medians move by at most 1.3 bootstrap standard errors; the 2100 95th percentile moves
    -3.6 cm (total) and -5.0 cm (Antarctica).
- **Land-water frame step** (`Ladrillo.LWS_OBS_ANCHOR = :zero_at_first_year`).
  - The observed land-water series is in the calibration target's 1995-2005 frame, while the other components start
    at 0 in 1850. v0.1 started the land-water component at the series' 1900 value, +1.565 cm, and Antarctica's
    sea-level feedback saw that step. v0.2 starts it at 0.
  - Glaciers, Greenland and thermal expansion are bit-identical.
  - Antarctica moves by about -1e-4 of its own contribution: at most 0.014 cm at 2100 and 0.049 cm at 2300 under
    SSP5-8.5.
- `--v0.1` on `scripts/project_joint.jl`, `project_fixed.jl` and `brick20_joint.jl` (or `posterior(paleo=:per_selection, ...)`
  and `LadrilloModel(...; lws_anchor=:obs_at_first_year)`) reproduces v0.1 bit for bit. The identity gate checks this.
- `test/identity/fix_effects.jl` separates the two fixes and checks each against what it should and should not move.

**v0.3** (2026-10-08). Two changes since `v0.2.0`, both made in the research code first and
certified bit-for-bit against it (CHANGELOG 2026-10-08i-k):
- **The van Vuuren forcing** (`data/forcing/*vv*`) is built from the PUBLISHED Zenodo v1.1.1 scenario. v0.2's files put
  the published post-2100 extension on the 2100 level of an earlier, gated prerelease; for Medium-to-Low that left the
  CO2 tail 10.5 GtCO2/yr low for two centuries. The SSP forcing is unchanged.
- **The joint-climate arm is conditioned on the record by default**: `scripts/project_joint.jl` drops every draw whose
  Antarctic fast dynamics fires by `record_end()` on its OWN paired FaIR config (3 of the 2,000 SSP draws, 2 of the
  van Vuuren). That is exact rejection, the same as the calibration likelihood's own weighting.
  `--no-record-conditioning` restores the v0.2 joint arm (outputs suffixed `_uncond`); `--v0.1` does not condition.
  The fixed-climate products are unchanged: on the mean climate no projection row crosses.
- `posterior(; drop_record_crossings=true)` drops the draws whose Antarctic fast dynamics fires by the end of the
  Antarctic record (`record_end()`, 2025) on the calibration's own climate. That is rejection sampling, i.e. exact
  conditioning of the joint (posterior x paleo) draws on "no fast dynamics through 2025". It drops 3 of the 10,000 draws
  (`record_last_year=2026`: 4), none of them among the 2,000 projection rows.
- `posterior(; drop_record_crossings=true)` tests ONE climate. In the joint-climate arm each draw runs on its own FaIR
  config, and more draws cross there; `scripts/project_joint.jl` therefore tests each pair on its own GMST with
  `fastdyn_onset_year(row, gmst, years)`, AFTER pairing (the arm pairs draws with configs by position).
- `fastdyn_onset_year` is the model's own trigger with the model's arithmetic. `test/conditioning/check_record_crossing.jl`
  checks it against model runs (10,000 draws on the record climate, 2,000 joint draws) and at the threshold to 1 ulp,
  and is mutation-tested. `scripts/record_crossings.jl` lists the onset year of every projection draw in each climate.

**v0.4** (2026-10-15). The BRICK 2.0 comparison arm reads the BRICK parameters of the
joint SNEASY-BRICK posterior published with MimiBRICK v2.0.0 (Zenodo 10.5281/zenodo.20592337), which
`tools/fetch_brick20_posterior.sh` downloads; v0.3 used an earlier subsample with no public URL (*BRICK 2.0 posterior*).
- Every Ladrillo product is bit-identical to v0.3. Only BRICK 2.0 numbers move: its SSP2-4.5 total is about 4 % lower
  at 2100 and 2300, and Table 5's profiled-AR(1) ΔAIC goes from 88.9 to 95.5.
- Certified bit-for-bit against the research re-run (identity gates 2-5; CHANGELOG 2026-10-15).

**v1.0** (the release the GMD paper reports). v0.4's model, posterior and numbers. Adds `python/figures/`, which
regenerates Figs. 1–6 byte-identical to the paper's renders, and `data/comparison/`, the FACTS and MAGICC-SLR
outputs those figures show. The options that reproduce v0.1 are renamed (`--v1.0` → `--v0.1`, `paleo=:v1` →
`paleo=:per_selection`, `lws_anchor=:v1_step` → `lws_anchor=:obs_at_first_year`); the old names stop with an error.
The internal tags `v1.0.0` to `v1.3.0` are renamed `v0.1.0` to `v0.4.0`.

## Known issues

- **The fast-dynamics threshold is crossed inside the hindcast by a few propagated draws.**
  - The calibration holds lambda and T_crit at fixed values that never cross (next item), but the paleo ensemble's
    low-T_crit tail does: 4 of the 10,000 draws (T_crit -16.4 to -16.9 C, DAIS scale) start fast dynamics in
    2021-2026 on the calibration climate, adding up to 1.6 cm of Antarctic sea level by 2026.
  - On its own FaIR config a draw can cross where it does not on the mean. In the joint-climate arm, 3 of the 2,000
    projection draws fire by 2025 under each SSP (2 under each van Vuuren marker), and 7-8 (5) by 2026; on the mean
    climate, none.
  - Hindcast products that subsample the posterior can therefore depend slightly on the paleo assignment.
  - v0.2 does not condition on the record. v0.3 conditions the joint-climate arm by default and offers the
    calibration-climate test as `posterior(; drop_record_crossings=true)` (*Versions*). The paper reports the
    conditioned numbers (v0.3 and later).
- **The calibration and the products use different values of the Antarctic gamma.**
  - The calibration holds lambda, T_crit and gamma at the central row's values (0.01365, -15.14 and 2.038;
    `Calibration.HELD`). The research calibrator was meant to hold them at their paleo medians, but it loaded the central
    row over them.
  - lambda and T_crit do not affect the hindcast. gamma does: `posterior()` attaches its paleo median, 2.834, so the
    projections and the posterior predictive use a gamma the chain was not conditioned on.
  - Measured on 500 draws: the hindcast moves by at most 0.015 cm; the projections by at most 0.13 cm through 2100, and
    by up to 3.2 cm in the 2300 Antarctic median (SSP5-8.5).
  - A calibration at gamma = 2.834 would give nearly the same posterior. Importance weights that convert the 10,000
    draws to it (the log-posterior difference between the two gammas) have an effective sample size of 99.7 %, and the
    reweighted SSP2-4.5 and SSP5-8.5 quantiles at 2100 and 2300 move by at most 0.54 cm. The products are therefore what
    a calibration at the paleo median would give, to that precision.

## License

Code: MIT (`LICENSE`).
