# data/comparison: static comparator outputs (FACTS, MAGICC-SLR)

These files hold the outputs of the two external models that Figs. 2, 3 and 6 compare Ladrillo and BRICK 2.0 against.
The package runs neither FACTS nor MAGICC-SLR; their outputs ship here as static data.
Each file is an unchanged copy of a file in the research repository SLR-RFF-BRICK, and `data/MANIFEST.csv` records its
source path, source commit and sha256. The package tests check every file against that sha256.

Readers in this package:
`python/figures/vv_model_comparison.py` builds `outputs/tables/vv_model_comparison.csv` from the first two files;
`python/figures/plot_model_comparison_components.py` (Figs. 2, 3) and `python/figures/plot_vv_responsiveness.py` (Fig. 6)
read that table, `comparator_classes.csv` and, for Fig. 6, `magicc_gmst_vv.csv`.

## Baselines

FACTS sea level is relative to its baseyear 2005. Every other sea level series in the comparison (Ladrillo, BRICK 2.0,
MAGICC-SLR) is relative to the 1995-2014 mean. FACTS sea level output has no years before 2005 (most modules start in
2020), so it cannot be re-referenced to 1995-2014; the comparison treats baseyear 2005 as approximately the 1995-2014
reference (the midpoint of that window is 2004.5), and the figure captions state the difference. Temperatures are relative to 1850-1900.

## Files

| File | What it is | Units | Made by (SLR-RFF-BRICK unless noted) |
|---|---|---|---|
| `facts_components_shared_n200.csv` | FACTS sea level quantiles (5/17/50/83/95 %) per scenario, component, module and year, n = 200 samples. Per-component modules (`ar5glaciers`, `emuglaciers`, `FittedISMIP`, `emuGrIS`, `ar5AIS`, `larmip`, `deconto21`, `emuAIS`, `bamber19`, `tlm`, `ssp-lws`) and workflow totals (`wf1f`, `wf2f`, `wf3f`, `wf4`, and the emulandice totals `wf1e`, `wf2e`, `wf3e`). The emulandice modules (`emu*`, `wf*e`) end at 2100; the others reach 2300. FACTS is driven by the same FaIR 2.2.4 (calib 1.6.0) CMIP7-basis GSAT as the Ladrillo and BRICK 2.0 projections; the `fair224-calib160` rows are that injected GSAT. Seven van Vuuren scenarios plus ssp126/245/585. | cm rel. FACTS baseyear 2005; GSAT rows degC rel. 1850-1900 | `extract_facts_shared_components.py` in the FACTS checkout, written to `outputs/facts_components_shared_n200.csv`; the 2300 rows come from separate 2300-horizon FACTS experiments merged under the same scenario names |
| `magicc_nauels_components_vv.csv` | MAGICC-SLR (MAGICC v7.5.3 with the Nauels 2025 sea level module) quantiles per scenario, component and year, 600-member AR6 drawnset. MAGICC computes its climate from the van Vuuren emissions, so this is the one comparator not on the shared FaIR driver. Column `arm`: `comparison` for the seven van Vuuren scenarios, `control` for ssp126/245/585, which re-run an earlier MAGICC result to check the pipeline. | cm rel. 1995-2014 | `python/extract_magicc_vv_components.py` |
| `magicc_gmst_vv.csv` | MAGICC's GMST quantiles per scenario and year from the same run, used in Fig. 6 for the High minus Very Low temperature gap on MAGICC's climate. | K rel. 1850-1900 | `python/extract_magicc_vv_gmst.py` |
| `comparator_classes.csv` | Classifies FACTS modules as `sej` (structured expert judgement: `bamber19`, and `wf4`, which uses bamber19 for both ice sheets) or `model` (every module not listed). The figures draw `sej` modules with an open marker. Lines starting with `#` are comments. | none | hand-maintained in the research repository as `benchmark/comparator_classes.csv` |

The FACTS run configurations, copies of the extraction scripts named above, and the receipt for the FACTS Greenland
numbers in Sect. 4.3 are in `comparison/` (see its README).

## Caption note

The Fig. 2 and 3 captions keep a sentence from the research renders saying the comparators are frozen for the SSP set
only. In this package every comparator file is pinned by its sha256 in `data/MANIFEST.csv`.
