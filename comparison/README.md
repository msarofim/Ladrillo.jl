# comparison: how the FACTS and MAGICC-SLR outputs were made

The paper compares Ladrillo with FACTS and MAGICC-SLR (Figs. 2, 3 and 6, and Sect. 4.3). This package does not run
either model; their outputs ship as static data in `data/comparison/`. The files here record how those outputs were
produced: the FACTS run configurations and run scripts, and the scripts that extracted both models' outputs.
They are not executed by this package and its tests.

Every file is an unchanged copy, taken from the commit that produced the data. `MANIFEST.csv` lists each file's source
repository, commit, path and sha256. `python tools/copy_comparison_scripts.py --check` verifies the copies against it.
Paths inside the scripts follow the authors' directory layout.

## FACTS (`facts/`)

- **Model.** FACTS v1.1.4 (https://github.com/radical-collaboration/facts, tag v1.1.4), with the container from
  `docker/Dockerfile`. The files here were added in a local checkout on top of that tag; commit `d52b9ed4`
  (2026-10-08) is the state that produced `data/comparison/facts_components_shared_n200.csv`.
- **Climate.** `build_shared_climate_nc.py` writes the FaIR 2.2.4 (calib 1.6.0) temperature ensemble that drives
  Ladrillo and BRICK 2.0 into FACTS's climate-input format, so all three models share one climate driver.
- **Configurations.** `build_shared_configs.py` generated `experiments/global.shared.<scenario>.n200/config.yml`
  (200 samples, to 2150) and the `<scenario>2300` experiments (to 2300), for the seven van Vuuren scenarios and
  SSP1-2.6, SSP2-4.5 and SSP5-8.5. The paper uses the van Vuuren scenarios.
- **Runs.** The `run_*.sh` scripts run the experiments in the container. The van Vuuren experiments were last run by
  `run_vv_zenodo_rerun_20261008.sh`.
- **Extraction.** `extract_facts_shared_components.py` reads the FACTS outputs and writes the quantile table. It
  appends each `*2300` experiment's years after 2150 to its base scenario, after checking that the two agree exactly on
  every overlapping year.
- **Greenland receipt (Sect. 4.3).** `experiments/global.shared.vvHL2300.n200.crate0` is the High-to-Low 2300
  experiment with one change, FittedISMIP's `crateyear_end: 0`, which evaluates the fit through 2300 instead of
  continuing each sample's 2080–2100 rate. This config was not committed in the FACTS checkout. `receipt/` holds the
  script that compares the two runs and its output: Greenland medians of 46 cm (default) and 112 cm (fit evaluated
  through 2300) at 2300, relative to 2005.

## MAGICC-SLR (`magicc/`)

- **Model.** MAGICC v7.5.3 with the sea level module of Nauels et al. (2025), on the 600-member AR6 parameter set with
  sea level parameters. It was run in the probabilistic MAGICC repository, which is available on request from the MAGICC
  team; its code and parameter set are not redistributed here.
- **Extraction.** `extract_magicc_vv_components.py` and `extract_magicc_vv_gmst.py` read that run's output and wrote
  `data/comparison/magicc_nauels_components_vv.csv` and `data/comparison/magicc_gmst_vv.csv`. Both import their
  conventions (baseline window, units, component map, quantiles) from `extract_magicc_components.py`.
