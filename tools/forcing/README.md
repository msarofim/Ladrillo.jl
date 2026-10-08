# tools/forcing — rebuild Ladrillo's FaIR forcing from its inputs

Regenerates all **43** files in `data/forcing/` — FaIR 2.2.4 (calib 1.6.0) GMST
(K) and ocean heat content (10²² J), anomalies vs 1850–1900, 841 configs — and
proves the rebuild matches the shipped files.

**Status (2026-10-08j, macOS arm64, Python 3.14.3, `requirements-forcing.txt`):**
all 43 files reproduce **byte-for-byte** (sha256 identical) from the inputs
below, with `--fairtofredi-root`, `--jobs 2`. The 28 van Vuuren (`vv*`) files
are built from the **published** Zenodo 20713982 v1.1.1 emissions; no gated
input is involved. All 23 intermediate CSVs the research repo also keeps
(calibration directories, emissions files, the 7 van Vuuren emissions files
included) match its copies byte-for-byte; the 7 van Vuuren splices match from
2023.5 on (their unused history columns differ by construction, see the chain).

**What a reviewer can rebuild from public inputs alone: 41 of the 43 files**
(measured 2026-10-08: an `--inputs-dir` holding every input except the one
non-public file, the AR6 SSP2-4.5 file; all 41 byte-identical). The two it
cannot are `fair_mean_{gmst,ohc}_ssp245harm.csv`, the calibration forcing,
because its emissions builder reads the AR6 file's post-2100 ratios.

| shipped file(s) | built by | scenario / forcing |
|---|---|---|
| `fair_mean_{gmst,ohc}_ssp245harm.csv` | `build_fair_mean_v160.py` | the **calibration** forcing: CMIP7 history + chrisroadmap SSP2-4.5, prescribed land use from CMIP7 marker M |
| `fair_cube_{gmst,ohc}_<ssp>_raw.csv`, `fair_mean_{gmst,ohc}_<ssp>.csv` for ssp126 / ssp245 / ssp585 | `build_fair_cube_v160.py` | CMIP7 history + RCMIP future; markers L / M / H |
| `fair_cube_gmst_ssp585_spliced.csv` | `build_fair_cube_spliced_v160.py` | mean path to 2014, then each config's anomaly re-referenced to 1995–2014 |
| `fair_cube_{gmst,ohc}_vv<M>_raw.csv`, `fair_mean_{gmst,ohc}_vv<M>.csv` for VL LN L ML M HL H | `build_fair_cube_vv_v160.py` | van Vuuren ScenarioMIP-CMIP7 markers: CMIP7 history 1750.5–2023.5 verbatim + the published Zenodo 20713982 v1.1.1 scenario 2024.5–2300.5, scaled per species at 2023.5; each marker with its own land use; `stochastic_run=False` |

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
#   optionally place the one NON-PUBLIC input (below) at the path it prints
python fetch_inputs.py --inputs-dir ~/forcing_inputs --verify     # sha256 vs inputs.sha256

python build_forcing.py --work-dir /tmp/forcing_build --inputs-dir ~/forcing_inputs --jobs 3
python check_forcing_identity.py --candidate-dir /tmp/forcing_build/forcing
#   without the non-public input: add --allow-missing (gates the 41 files built)
```

In the maintainer's research checkout every input is already present, so the
fetch can be skipped:

```bash
python build_forcing.py --work-dir /tmp/forcing_build \
    --fairtofredi-root <research-checkout> --jobs 3
```

`build_forcing.py` writes only under `--work-dir`. It touches `data/forcing/`
**only** with `--install`, and only once all 43 files were built. Before it runs
anything it checks each input against `inputs.sha256` (RCMIP also against fair's
own md5, in pooch's cache, so fair cannot download it mid-build):

- an input that is **present but fails** its checksum **aborts the build**;
  nothing runs;
- an input that is **absent** drops the steps that read it and every step
  downstream of them. The build prints which steps and which of the 43 files it
  will not build, then builds the rest; `--install` refuses.

Each step declares the inputs it reads (the calibration setup steps take a
directory, so their files come from the scripts' own file-name constants, and
step 4 reads RCMIP through `fair.fill_from_rcmip`). It runs the builders from a
frozen copy (`<work-dir>/tools_snapshot/`), pins BLAS to one thread, and writes
`<work-dir>/BUILD_PROVENANCE.json`: every command and its wall time, the sha256 of
each script and input, the inputs that were absent and the steps and files that
were therefore skipped, package versions, and the Ladrillo commit.

### Runtime

Measured on a 10-core Apple-silicon Mac:

| stage | load average about 3 | load average about 35 (2026-10-08j) |
|---|---|---|
| calib + emissions (6 steps, no FaIR) | ~9 s | 16 s |
| each FaIR run (841 configs, 1750–2301, single-threaded) | 12–17 s | 22–33 s |
| **full build** | **55–69 s** (`--jobs 3`, two runs) | **179 s** (`--jobs 2`); 147 s without the AR6 file |
| full build, `--jobs 1` | 1 min 56 s | — |
| identity gate (43 files) | < 1 s | — |

Each FaIR run needs roughly 1–2 GB of RAM.

## The chain

| # | builder | reads | writes (under `--work-dir`) |
|---|---|---|---|
| 1 | `setup_calibration_v160.py` | Zenodo 18828694 params, species, volcanic, solar, irrigation | `calibration_v160/` (Land use → `calculated`, shared irrigation; this is the VV-vs-RFF comparison stack) |
| 2 | `setup_calibration_v160_production.py` | 1 + Zenodo land use, volcanic, solar, irrigation | `calibration_v160_prod/` (stock 1.6.0 land use; `volcanic_solar_<MARKER>.csv` ×7) |
| 3 | `build_emissions_v160_ssp245.py` | CMIP7 history, chrisroadmap SSP2-4.5, AR6 r(t) | `calibration_v160/emissions_v160_cmip7_ssp245_harmonized.csv` |
| 4 | `build_emissions_v160_ssp_rcmip.py` | CMIP7 history, RCMIP v5.1.0 (via `fill_from_rcmip`) | `calibration_v160_prod/emissions_v160_cmip7harm_{ssp126,ssp245,ssp585}.csv` |
| 6 | `build_vanvuuren_extended_splice.py --skeleton` | Zenodo 20713982 v1.1.1, `vv_template_skeleton.csv` | `vanvuuren/spliced_ext_zenodo/` |
| 7 | `build_emissions_v160_cmip7harm_vv.py` | 6 + CMIP7 history | `calibration_v160_prod/emissions_v160_cmip7harm_vv<M>.csv` |
| 8 | `build_fair_mean_v160.py` | 2 + 3 | `forcing/fair_mean_*_ssp245harm.csv` |
| 9 | `build_fair_cube_v160.py` ×3 | 2 + 4 | `forcing/fair_{cube,mean}_*_<ssp>*` |
| 10 | `build_fair_cube_vv_v160.py` ×7 | 2 + 7 | `forcing/fair_{cube,mean}_*_vv<M>*` |
| 11 | `build_fair_cube_spliced_v160.py` | 9 (ssp585) | `forcing/fair_cube_gmst_ssp585_spliced.csv` |

**There is no step 5 any more** (removed 2026-10-08j; the step ids keep their
numbers). It was `build_vanvuuren_spliced.py`, which built a van Vuuren template
from the gated IIASA `scenariomip-2026-prerelease` file plus Smith 2024 history,
and step 6 filled that template. On the `zenodo` variant every value from 2023.5
on comes from Zenodo, and step 7 reads only the 2023.5 anchor and the later
columns (history is the CMIP7 record). The template therefore contributed only
row order, the `scenario` / `region` / `variable` / `unit` strings and the
column grid. `vv_template_skeleton.csv` carries exactly those names and units
(no numbers), taken from the research repo's `data/vanvuuren/spliced_ext_zenodo/`
files, one row per marker × species, in the template's row order. With
`--skeleton`, step 6 reproduces the template's grid, including its gap (2021.5
and 2022.5 are absent), and fills **every** column from Zenodo. Its history
columns therefore carry Zenodo's own history instead of Smith 2024; nothing
downstream reads them. The step 7 outputs are byte-identical either way.
`build_vanvuuren_spliced.py` stays in this directory for the template-based
variants (`--spliced-dir`), which still need the IIASA file.

Helpers: `fair_basis.py`, `fair_units.py` (verbatim copies), `forcing_inputs.py`
(the input table), and two vendored constant modules, `v145_species_maps.py` and
`emissions_v145_utils.py`. Those two are copied from private research-repo
scripts as source text, and checked equal to the originals.

## Inputs

`inputs.sha256` lists all 13 files the builders read, by their path under
`--inputs-dir`. The table below says where each comes from, under what licence,
and what it feeds. Licences were read from the Zenodo and GitHub APIs (metadata
only) on 2026-10-08.

| input (`--inputs-dir` path) | public? | source | licence | feeds |
|---|---|---|---|---|
| `zenodo_18828694/calibrated_constrained_parameters.csv` | yes | Zenodo [10.5281/zenodo.18828694](https://doi.org/10.5281/zenodo.18828694), fair-calibrate v1.6.0 (2026-03-02). The local copy is byte-identical to a pooch download from `zenodo.org/records/18828694/files/…` | CC-BY-4.0 (attribution below) | every file |
| `zenodo_18828694/species_configs_properties.csv` | yes | same record | CC-BY-4.0 | every file |
| `zenodo_18828694/volcanic_forcing_timebounds_cmip7.csv` | yes | same record | CC-BY-4.0 | every file |
| `zenodo_18828694/solar_forcing_timebounds_cmip7.csv` | yes | same record | CC-BY-4.0 | every file |
| `zenodo_18828694/irrigation_forcing_timebounds_cmip7.csv` | yes | same record | CC-BY-4.0 | every file |
| `zenodo_18828694/land_use_forcing_timebounds_cmip7.csv` | yes | same record | CC-BY-4.0 | every file |
| `zenodo_18828694/historical_emissions_1750-2023_cmip7.csv` | yes | same record | CC-BY-4.0 | every file |
| `rcmip/rcmip-emissions-annual-means-v5-1-0.csv` | yes | RCMIP v5.1.0, fetched by fair 2.2.4 through pooch (md5 pinned by fair) | not checked | ssp126/245/585 |
| `rcmip/rcmip-concentrations-annual-means-v5-1-0.csv` | yes | as above. fair reads it, but only emissions-driven species are defined, so its values are unused | not checked | — |
| `rcmip/rcmip-radiative-forcing-annual-means-v5-1-0.csv` | yes | as above, likewise unused | not checked | — |
| `chrisroadmap_fair-calibrate/CH4Pledge-rcmip-emissions-ssp245-CLE-MFR-filled.csv` | yes | GitHub chrisroadmap/fair-calibrate, `data/emissions/`, URL pinned to commit `3f72545` (the last commit touching the file, from GitHub API metadata; nothing was downloaded). ⚠ Whether that commit's bytes are the copy the shipped forcing was built from is **not yet verified**: `fetch_inputs.py --download` then `--verify` checks them against `inputs.sha256` on first use | MIT | ssp245harm, 2024–2100 |
| `zenodo_20713982/emissions_1750-2500.csv` | yes | Zenodo [10.5281/zenodo.20713982](https://doi.org/10.5281/zenodo.20713982) **v1.1.1** (2026-06-16), member `data/fair-inputs/emissions_1750-2500.csv` of `chrisroadmap/cmip7-scenariomip-1.1.1.zip`. ⚠ Not the DOI printed in van Vuuren et al. 2026, which resolves to the stylized v1.0 | MIT | **all 28 vv files**: every value from 2023.5 on (the anchor and 2024.5–2300.5) |
| `edf_rennels/AR6_emissions_ssp245_1750_2300.csv` | ⛔ **NO** | EDF / Lisa Rennels' AR6 SSP2-4.5 file, as received. No public origin is recorded | none recorded | ssp245harm, **2101–2300 only**: its ratios r(t) shape the tail; a perturbation moved no emissions column ≤ 2100.5 |

**Attribution (CC-BY-4.0, Zenodo 18828694).** Smith, C. (2026). *fair
calibration data* (v1.6.0) [Data set]. Zenodo.
https://doi.org/10.5281/zenodo.18828694. The calibration directories this build
writes, and every forcing file in `data/forcing/`, are derived from it.

### Flags

- **Not public: one input**, the Rennels AR6 file. Without it, a reviewer
  rebuilds 41 of the 43 files (measured, above) and cannot rebuild
  `ssp245harm`. The build cannot run step 3 at all without the file, although
  the file shapes only the 2101–2300 tail. The calibration window is 1850–2026,
  entirely before that tail; this follows from causality and was not measured on
  the FaIR output.
- **No longer inputs (2026-10-08j):** the IIASA `scenariomip-2026-prerelease`
  pull and the Smith 2024 history (chrisroadmap/temperature-attribution). Both
  fed only the van Vuuren template that `vv_template_skeleton.csv` replaces.
- **No licence recorded:** the AR6 file. **Not checked:** RCMIP.
- **Derived locally, now regenerated:** `calibration_v160/` and
  `calibration_v160_prod/` (untracked in the research repo, no LICENCE), every
  `emissions_v160_*` file, and `vanvuuren/spliced_ext_zenodo/`. They are not
  inputs. Steps 1–4, 6 and 7 rebuild them from the files above.
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

There is **no tolerance**. A missing file FAILS. `--allow-missing` gates a
partial build instead: every file that exists must pass as above, the absent
ones are listed as ABSENT and do not fail, and at least one must exist.
`--report FILE` writes the table as CSV. `--self-test` mutation-tests the gate
itself. It must FAIL a one-ulp change to one cell (it does: 2.2e-16) and a
dropped row. It must PASS on values, and FAIL on sha256, a CRLF-only change.

**Mutation tests.** Rebuilding the calibration forcing with the plausible wrong
recipe, `--calib-dir calibration_v160 --marker none` (2026-10-08), gives:

| file | sha256 | values | max \|diff\| | verdict |
|---|---|---|---|---|
| `fair_mean_gmst_ssp245harm.csv` | no | no | 0.0372 K (2055); −0.0302 K at 2024 | FAIL |
| `fair_mean_ohc_ssp245harm.csv` | no | no | 1.62 ×10²² J (1959) | FAIL |

A second mutation, `vvM` with `--stochastic on` (the SSP-cube convention), FAILS
all four vvM files: 0.65 K max per config, 0.015 K on the mean. 2026-10-08j: one
digit changed in a built file FAILS the gate (`fair_cube_ohc_vvML_raw.csv`, max
\|diff\| 5.0e-4, with or without `--allow-missing`), and one digit
changed in a copy of an input aborts `build_forcing.py` before any step runs,
with or without another input absent.

`--intermediates WORK --fairtofredi-root ROOT` compares a build's intermediates
with the research repo's copies (maintainer's research checkout only). For
`vanvuuren/spliced_ext_zenodo/` it compares the names and the columns from
2023.5 on; the history columns differ by construction (see the chain).

## What the port changed

These are ports, not rewrites. Each builder's numerical code is copied from
FaIRtoFrEDI @ `8984a13`, and the only changes are to I/O:

- Hard-coded `PROJ_ROOT` / research-checkout paths became required CLI flags.
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
  here reads, is dropped. Since 2026-10-08j the chain does not run it.
- `build_vanvuuren_extended_splice.py` gained `--skeleton` (2026-10-08j), a
  Ladrillo addition: the `zenodo` variant without the template.
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
