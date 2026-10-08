# tools/forcing — rebuild Ladrillo's FaIR forcing from its inputs

Regenerates all **43** files in `data/forcing/` — FaIR 2.2.4 (calib 1.6.0) GMST
(K) and ocean heat content (10²² J), anomalies vs 1850–1900, 841 configs — and
proves the rebuild matches the shipped files.

**Status (2026-10-08, macOS arm64, Python 3.14.3, `requirements-forcing.txt`):**
all 43 files reproduce **byte-for-byte** (sha256 identical), on three rebuilds
(`--fairtofredi-root` with the research-repo gates on; `--inputs-dir` with
`--jobs 3`; `--inputs-dir` with `--jobs 1`). Every intermediate (calibration directories,
emissions files, van Vuuren splices) also matches the research repo's copy
byte-for-byte.

| shipped file(s) | built by | scenario / forcing |
|---|---|---|
| `fair_mean_{gmst,ohc}_ssp245harm.csv` | `build_fair_mean_v160.py` | the **calibration** forcing: CMIP7 history + chrisroadmap SSP2-4.5, prescribed land use from CMIP7 marker M |
| `fair_cube_{gmst,ohc}_<ssp>_raw.csv`, `fair_mean_{gmst,ohc}_<ssp>.csv` for ssp126 / ssp245 / ssp585 | `build_fair_cube_v160.py` | CMIP7 history + RCMIP future; markers L / M / H |
| `fair_cube_gmst_ssp585_spliced.csv` | `build_fair_cube_spliced_v160.py` | mean path to 2014, then each config's anomaly re-referenced to 1995–2014 |
| `fair_cube_{gmst,ohc}_vv<M>_raw.csv`, `fair_mean_{gmst,ohc}_vv<M>.csv` for VL LN L ML M HL H | `build_fair_cube_vv_v160.py` | van Vuuren ScenarioMIP-CMIP7 markers on CMIP7 history, each with its own marker's land use; `stochastic_run=False` |

## ⚠ The calibration forcing is a MIX of two directories

`ssp245harm` is `calibration_v160_prod` (land use **prescribed**, marker **M**)
driven on `calibration_v160`'s emissions file. The obvious
`--calib-dir calibration_v160 --marker none` lands **0.030 K cooler at 2024**
(max 0.037 K, in 2055) and reproduces neither file — that is the mutation test
below. The emissions file alone does not determine the output.
`build_forcing.py` encodes the right combination; don't hand-assemble it.

## Quick start

```bash
python -m venv .venv-forcing && . .venv-forcing/bin/activate
pip install -r tools/forcing/requirements-forcing.txt
cd tools/forcing

python fetch_inputs.py --inputs-dir ~/forcing_inputs              # plan only, no network
python fetch_inputs.py --inputs-dir ~/forcing_inputs --download   # fetch the PUBLIC inputs
#   then place the two NON-PUBLIC inputs (below) at the paths it prints
python fetch_inputs.py --inputs-dir ~/forcing_inputs --verify     # sha256 vs inputs.sha256

python build_forcing.py --work-dir /tmp/forcing_build --inputs-dir ~/forcing_inputs --jobs 3
python check_forcing_identity.py --candidate-dir /tmp/forcing_build/forcing
```

On Marcus's machine the research repo already holds every input, so skip the fetch:

```bash
python build_forcing.py --work-dir /tmp/forcing_build \
    --fairtofredi-root ~/Documents/2026/CodeProjects/FaIRtoFrEDI --jobs 3
```

`build_forcing.py` writes only under `--work-dir`. It touches `data/forcing/`
**only** with `--install`, and only once all 43 files exist. Before it runs
anything it checks every input against `inputs.sha256` (and that RCMIP is already
in pooch's cache, so fair cannot download it mid-build). It runs the builders from
a frozen copy (`<work-dir>/tools_snapshot/`), pins BLAS to one thread, and writes
`<work-dir>/BUILD_PROVENANCE.json`: every command and its wall time, the sha256 of
each script and input, package versions, and the Ladrillo commit.

### Runtime

Measured on a 10-core Apple-silicon Mac, load average about 3:

| stage | wall |
|---|---|
| calib + emissions (7 steps, no FaIR) | ~9 s |
| each FaIR run (841 configs, 1750–2301, single-threaded) | 12–17 s |
| **full build, `--jobs 3`** | **55–69 s** (two runs) |
| full build, `--jobs 1` | 1 min 56 s |
| identity gate (43 files) | < 1 s |

Each FaIR run needs roughly 1–2 GB of RAM.

## The chain

| # | builder | reads | writes (under `--work-dir`) |
|---|---|---|---|
| 1 | `setup_calibration_v160.py` | Zenodo 18828694 params, species, volcanic, solar, irrigation | `calibration_v160/` (Land use → `calculated`, shared irrigation; this is the VV-vs-RFF comparison stack) |
| 2 | `setup_calibration_v160_production.py` | 1 + Zenodo land use, volcanic, solar, irrigation | `calibration_v160_prod/` (stock 1.6.0 land use; `volcanic_solar_<MARKER>.csv` ×7) |
| 3 | `build_emissions_v160_ssp245.py` | CMIP7 history, chrisroadmap SSP2-4.5, AR6 r(t) | `calibration_v160/emissions_v160_cmip7_ssp245_harmonized.csv` |
| 4 | `build_emissions_v160_ssp_rcmip.py` | CMIP7 history, RCMIP v5.1.0 (via `fill_from_rcmip`) | `calibration_v160_prod/emissions_v160_cmip7harm_{ssp126,ssp245,ssp585}.csv` |
| 5 | `build_vanvuuren_spliced.py` | IIASA prerelease 2023–2100, Smith 2024 | `vanvuuren/spliced/` |
| 6 | `build_vanvuuren_extended_splice.py` | 5 + Zenodo 20713982 v1.1.1 | `vanvuuren/spliced_ext_harmonized/` |
| 7 | `build_emissions_v160_cmip7harm_vv.py` | 6 + CMIP7 history | `calibration_v160_prod/emissions_v160_cmip7harm_vv<M>.csv` |
| 8 | `build_fair_mean_v160.py` | 2 + 3 | `forcing/fair_mean_*_ssp245harm.csv` |
| 9 | `build_fair_cube_v160.py` ×3 | 2 + 4 | `forcing/fair_{cube,mean}_*_<ssp>*` |
| 10 | `build_fair_cube_vv_v160.py` ×7 | 2 + 7 | `forcing/fair_{cube,mean}_*_vv<M>*` |
| 11 | `build_fair_cube_spliced_v160.py` | 9 (ssp585) | `forcing/fair_cube_gmst_ssp585_spliced.csv` |

Helpers: `fair_basis.py`, `fair_units.py` (verbatim copies), `forcing_inputs.py`
(the input table), and two vendored constant modules, `v145_species_maps.py` and
`emissions_v145_utils.py`. Those two are copied from private research-repo
scripts as source text, and checked equal to the originals.

## Inputs

`inputs.sha256` lists all 15 files the builders read, by their path under
`--inputs-dir`. The table below says where each comes from and what it feeds.

| input (`--inputs-dir` path) | public? | source | feeds |
|---|---|---|---|
| `zenodo_18828694/calibrated_constrained_parameters.csv` | yes | Zenodo [10.5281/zenodo.18828694](https://doi.org/10.5281/zenodo.18828694), fair-calibrate v1.6.0 (2026-03-02). The local copy is byte-identical to a pooch download from `zenodo.org/records/18828694/files/…` | every file |
| `zenodo_18828694/species_configs_properties.csv` | yes | same record | every file |
| `zenodo_18828694/volcanic_forcing_timebounds_cmip7.csv` | yes | same record | every file |
| `zenodo_18828694/solar_forcing_timebounds_cmip7.csv` | yes | same record | every file |
| `zenodo_18828694/irrigation_forcing_timebounds_cmip7.csv` | yes | same record | every file |
| `zenodo_18828694/land_use_forcing_timebounds_cmip7.csv` | yes | same record | every file |
| `zenodo_18828694/historical_emissions_1750-2023_cmip7.csv` | yes | same record | every file |
| `rcmip/rcmip-emissions-annual-means-v5-1-0.csv` | yes | RCMIP v5.1.0, fetched by fair 2.2.4 through pooch (md5 pinned by fair) | ssp126/245/585 |
| `rcmip/rcmip-concentrations-annual-means-v5-1-0.csv` | yes | as above. fair reads it, but only emissions-driven species are defined, so its values are unused | — |
| `rcmip/rcmip-radiative-forcing-annual-means-v5-1-0.csv` | yes | as above, likewise unused | — |
| `chrisroadmap_fair-calibrate/CH4Pledge-rcmip-emissions-ssp245-CLE-MFR-filled.csv` | yes, **URL unverified** | GitHub chrisroadmap/fair-calibrate `data/emissions/`. The path comes from research notes (2026-05-27). The branch in the URL, and whether the fetched bytes match, are **not verified** (no download was made); the sha256 check will decide | ssp245harm, 2024–2100 |
| `chrisroadmap_temperature-attribution/emissions_extrapolated_1750-2025.csv` | yes | GitHub chrisroadmap/temperature-attribution @ `4b6022a`, `output/` (MIT). sha256 verified against a local clone | vv files: **species list and unit labels only**. A perturbation of every Smith value left all 7 vv emissions files byte-identical |
| `zenodo_20713982/emissions_1750-2500.csv` | yes | Zenodo [10.5281/zenodo.20713982](https://doi.org/10.5281/zenodo.20713982) **v1.1.1** (2026-06-16), member `data/fair-inputs/emissions_1750-2500.csv` of `chrisroadmap/cmip7-scenariomip-1.1.1.zip`. ⚠ Not the DOI printed in van Vuuren et al. 2026, which resolves to the stylized v1.0 | vv files, 2101–2300 |
| `iiasa_scenariomip-2026-prerelease/emissions_harmonized_infilled_2023-2100.csv` | ⛔ **NO** | IIASA Scenario Explorer `scenariomip-2026-prerelease`, harmonized + infilled, 7 markers. **Gated (login) prerelease**; redistribution terms not established | **all 28 vv files**, 2024–2100 (perturbing it changes all 7 vv emissions files) |
| `edf_rennels/AR6_emissions_ssp245_1750_2300.csv` | ⛔ **NO** | EDF / Lisa Rennels' AR6 SSP2-4.5 file, as received. No public origin or licence is recorded | ssp245harm, **2101–2300 only**: its ratios r(t) shape the tail; a perturbation moved no emissions column ≤ 2100.5 |

### Flags

- **Not public: two inputs.** Without the IIASA prerelease file, a reviewer cannot
  rebuild any `vv*` file. Without the Rennels AR6 file, a reviewer cannot rebuild
  `ssp245harm` past 2100. The calibration window is 1850–2026, which lies entirely
  before that tail; this follows from causality and was not measured on the FaIR
  output. Shipping either file, or a public substitute, is a decision for Marcus.
  The build does not reconstruct either one. Without them, its input check stops
  it.
- **No licence recorded:** the AR6 file. The licences of the two Zenodo records,
  RCMIP and fair-calibrate were not checked offline.
- **Derived locally, now regenerated:** `calibration_v160/` and
  `calibration_v160_prod/` (untracked in the research repo, no LICENCE), every
  `emissions_v160_*` file, and `data/vanvuuren/spliced*/`. They are no longer
  inputs. Steps 1–7 rebuild them byte-identically from the files above.
- **Unverified fetches.** Nothing was downloaded while this port was built. Only
  one Zenodo file is confirmed byte-identical to its public URL (the parameters,
  via pooch's cache name). For every other public input, `fetch_inputs.py
  --download` followed by `--verify` is the first real test that the public copy
  is the one used.

## The identity gate

`check_forcing_identity.py --candidate-dir DIR` lists the shipped forcing files
from `data/MANIFEST.csv`. It first checks each reference file against the
manifest's sha256, so a reference edited in place cannot pass. Then it compares
each candidate file:

1. sha256 equal → PASS;
2. otherwise both files are read with `float_precision="round_trip"`. The file
   passes only if column names and shape are identical and **every cell is `==`**
   (formatting-only differences pass here and are labelled as such);
3. anything else → FAIL, with the max |difference| and where it occurs.

There is **no tolerance**. `--report FILE` writes the table as CSV. `--self-test`
mutation-tests the gate itself. It must FAIL a one-ulp change to one cell (it does:
2.2e-16) and a dropped row. It must PASS on values, and FAIL on sha256, a
CRLF-only change.

**Mutation test (2026-10-08).** Rebuilding the calibration forcing with the
plausible wrong recipe, `--calib-dir calibration_v160 --marker none`, gives:

| file | sha256 | values | max \|diff\| | verdict |
|---|---|---|---|---|
| `fair_mean_gmst_ssp245harm.csv` | no | no | 0.0372 K (2055); −0.0302 K at 2024 | FAIL |
| `fair_mean_ohc_ssp245harm.csv` | no | no | 1.62 ×10²² J (1959) | FAIL |

A second mutation, `vvM` with `--stochastic on` (the SSP-cube convention), FAILS
all four vvM files: 0.65 K max per config, 0.015 K on the mean.

`--intermediates WORK --fairtofredi-root ROOT` compares a build's intermediates
with the research repo's copies (Marcus's machine only).

## What the port changed

These are ports, not rewrites. Each builder's numerical code is copied from
FaIRtoFrEDI @ `8984a13`, and the only changes are to I/O:

- Hard-coded `PROJ_ROOT` / `~/Documents/...` paths became required CLI flags.
- Every output defaults to a directory the caller names. None defaults to
  `data/forcing/`.
- `setup_calibration_v160.py` reads the Zenodo files from `--raw-dir` instead of
  downloading them. The download moved to `fetch_inputs.py`. The bytes still go
  through `io.BytesIO` as before.
- Research-repo-only **gates** became optional flags that print when they are
  skipped. They are the cube-builder `[CUBE-TRUNCATION]` check, the spliced-cube
  PORT GATE (`--port-gate-dir`; without it `--skip-port-gate` must be explicit),
  and the vendored-map cross-check (`--v145-builder`). All three PASS when they
  are pointed at the research repo.
- `build_vanvuuren_spliced.py` is lines 67–232 of the research
  `run_fair_vanvuuren_erf.py`, verbatim. Its v1.4.5 FaIR ERF run, which nothing
  here reads, is dropped.
- The SSP2-4.5 builder's closing printout used to read a v1.4.5 research file. It
  is now opt-in (`--compare-v145-file`) and never touched the output.

Each script's docstring has a `PORT` paragraph listing its own changes.

## Conventions worth knowing (not changed here)

- **The vv cubes are deterministic (`stochastic_run=False`); the SSP cubes keep
  the posterior's `stochastic_run=True`.** The difference is deliberate in the
  research builders (see `build_fair_cube_vv_v160.py`'s docstring) and the
  shipped files follow it. The second mutation above shows the gate would catch a
  mix-up.
- The SSP → CMIP7-marker mapping for prescribed land use and irrigation is
  ssp126→L, ssp245→M, ssp585→H. It is bit-neutral through 2023; the research repo
  measured 0.02–0.09 K at 2100/2300.
- Byte identity was shown on one platform and one package set. Another numpy or
  BLAS build may differ in the last bits. If it does, the gate reports exactly how
  much.
