"""
build_vanvuuren_spliced.py

The SPLICE half of the research repo's FaIRtoFrEDI/scripts/run_fair_vanvuuren_erf.py
@ 8984a13: it writes the per-marker wide-format van Vuuren emissions CSVs that
the research repo keeps in data/vanvuuren/spliced/<scenario_slug>.csv.

Splice convention (unchanged):
  1750.5 - 2020.5  : Smith 2024 historical (chrisroadmap/temperature-attribution)
                     [the research driver says "- 2022.5", but load_smith() stops
                     at SMITH_LAST_YEAR_INT = 2020, so 2021.5 and 2022.5 are
                     ABSENT from the grid -- the gap fair interpolates across,
                     noted in calibration_v160_prod/README.md]
  2023.5 - 2100.5  : van Vuuren ScenarioMIP-CMIP7 (IIASA scenariomip-2026-prerelease,
                     harmonized + infilled)
  2101.5 - 2300.5  : hold constant at each species' 2100.5 value (the placeholder
                     tail that build_vanvuuren_extended_splice.py replaces)

PORT (Ladrillo tools/forcing/, 2026-10-08). Lines 67-232 of the research driver
(its unit/species maps, slugify, Step 1 and Step 2) are copied VERBATIM; Steps
3-8 (the v1.4.5 FaIR ERF run and its .npz) are dropped because nothing in
Ladrillo's forcing reads them. Only the I/O plumbing differs: --vv-file,
--smith-file and --out-dir replace the research repo's hard-coded paths, and
`emissions_v145_utils` is the vendored subset in this directory.

⚠ THE INPUT IS NOT PUBLIC. emissions_harmonized_infilled_2023-2100.csv is a pull
from the GATED IIASA scenariomip-2026-prerelease database.

NO LONGER IN THE BUILD (2026-10-08). build_forcing.py does not run this script:
the shipped vv files are on the `zenodo` variant, which
build_vanvuuren_extended_splice.py builds from the published Zenodo file and
vv_template_skeleton.csv alone, so no IIASA value reaches a shipped file. This
script is kept for the template-based variants (`harmonized`, `native`,
`fullnative`; `--spliced-dir`), which do need the IIASA file. Neither of its two
inputs is listed in forcing_inputs.py / inputs.sha256 any more; pass them by path.
"""

from pathlib import Path
import argparse
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import emissions_v145_utils as _eu  # noqa: E402
from emissions_v145_utils import (  # noqa: E402
    SPECIES_MAP as RFF_SPECIES_MAP,
    RFF_TO_SMITH_SCALE,
    load_smith,
)


# Per-FaIR-species VV→Smith unit conversion. Merges RFF_TO_SMITH_SCALE
# (handles CO2 FFI/AFOLU Mt→Gt) with the VV-specific N2O kt→Mt fix. The
# 2026-prerelease reports N2O in `kt N2O/yr` while Smith 2024 uses
# `Mt N2O/yr`; without this the driver silently feeds N2O emissions 1000×
# too large, blowing up N2O concentrations and cascading through the CH4
# lifetime chemistry (log-of-negative → NaN for both CH4 and N2O ERF).
VV_TO_SMITH_SCALE = {
    **RFF_TO_SMITH_SCALE,
    "N2O": 1e-3,       # kt N2O/yr (VV) → Mt N2O/yr (Smith)
}

# ── Configuration ──────────────────────────────────────────────────────────────
START_YEAR = 1750
END_YEAR   = 2301
VV_LAST_YEAR = 2100          # last VV year (used at timepoint 2100.5)
SMITH_LAST_YEAR = 2022       # last Smith year kept (timepoint 2022.5)
TAIL_END_YEAR   = 2300       # hold constant to 2300.5

# PORT: the input file, the output directory and the Smith file come from the
# CLI. The v1.4.5 calibration paths and the ERF-cube output of the research
# driver are not needed for the splice and are dropped with the FaIR run.
_ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
_ap.add_argument("--vv-file", type=Path, required=True,
                 help="emissions_harmonized_infilled_2023-2100.csv (IIASA "
                      "scenariomip-2026-prerelease pull; NOT public)")
_ap.add_argument("--smith-file", type=Path, required=True,
                 help="Smith 2024 history: chrisroadmap/temperature-attribution "
                      "output/emissions_extrapolated_1750-2025.csv")
_ap.add_argument("--out-dir", type=Path, required=True,
                 help="written as <out-dir>/<scenario_slug>.csv (the research "
                      "repo's data/vanvuuren/spliced/)")
_args = _ap.parse_args()
_eu.SMITH_FILE = _args.smith_file

VV_FILE      = _args.vv_file
SPLICED_DIR  = _args.out_dir
SPLICED_DIR.mkdir(parents=True, exist_ok=True)

VV_PREFIX = "Climate Assessment|Harmonized and Infilled|Emissions|"

# Species dropped from VV (not represented in v1.4.5 or would double-count)
VV_DROP_SPECIES = {
    "Halon1202",              # v1.4.5 dropped it
    "CO2",                    # we use FFI + AFOLU split
    "GHG AR6GWP100",          # GWP aggregate
    "Kyoto GHG AR6GWP100",    # GWP aggregate
}

# Extra VV-only species mappings: VV files PFCs under bare names (CF4, C2F6, ...)
# while RFF-SP files them as PFC|CF4, PFC|C2F6, .... The bare names collide with
# the RFF splice's SPECIES_MAP keys and would silently drop if we used the
# RFF map alone. Extend it here.
VV_EXTRA_MAP = {
    "CF4":   "CF4",
    "C2F6":  "C2F6",
    "C3F8":  "C3F8",
    "cC4F8": "c-C4F8",
    "C4F10": "C4F10",
    "C5F12": "C5F12",
    "C6F14": "C6F14",
    "C7F16": "C7F16",
    "C8F18": "C8F18",
}
VV_SPECIES_MAP = {**RFF_SPECIES_MAP, **VV_EXTRA_MAP}


def slugify(name):
    s = name.lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_")


# ── Step 1: Load Smith historical + build the year column layout ──────────────
print(f"── van Vuuren splice (Smith 2024 history + IIASA 2023-2100, constant tail) ──")
print(f"Loading Smith 2024 historical ...")
df_smith_full, keep_smith_all = load_smith()
# We use Smith through 2022.5 (the year before VV starts, VV starts 2023).
keep_smith = [c for c in keep_smith_all if float(c) <= SMITH_LAST_YEAR + 0.5]
print(f"  Smith kept: {keep_smith[0]} ... {keep_smith[-1]}  ({len(keep_smith)} timepoints)")

vv_year_cols  = [f"{y}.5" for y in range(SMITH_LAST_YEAR + 1, VV_LAST_YEAR + 1)]
tail_year_cols = [f"{y}.5" for y in range(VV_LAST_YEAR + 1, TAIL_END_YEAR + 1)]
all_year_cols = keep_smith + vv_year_cols + tail_year_cols
print(f"  VV block:  {vv_year_cols[0]} ... {vv_year_cols[-1]}  ({len(vv_year_cols)})")
print(f"  Tail (constant): {tail_year_cols[0]} ... {tail_year_cols[-1]}  ({len(tail_year_cols)})")
print(f"  Total timepoints: {len(all_year_cols)}")

# ── Step 2: Load VV long-format CSV -> wide by species per scenario ──────────
print(f"\nLoading VV long-format from {VV_FILE.name} ...")
df_vv = pd.read_csv(VV_FILE)
scenarios = sorted(df_vv.scenario.unique().tolist())
print(f"  Scenarios ({len(scenarios)}):")
for s in scenarios:
    print(f"    - {s!r}")

# Strip prefix and map to canonical FaIR names.
df_vv["specie_suffix"] = df_vv["variable"].str.removeprefix(VV_PREFIX)
# Filter drops first (aggregates and Halon1202)
df_vv = df_vv[~df_vv["specie_suffix"].isin(VV_DROP_SPECIES)].copy()
df_vv["fair_specie"] = df_vv["specie_suffix"].map(VV_SPECIES_MAP)

unmapped = sorted(df_vv[df_vv["fair_specie"].isna()]["specie_suffix"].unique().tolist())
if unmapped:
    print(f"  WARN: {len(unmapped)} VV species not in SPECIES_MAP (dropped):")
    for u in unmapped:
        print(f"    - {u}")
df_vv = df_vv[df_vv["fair_specie"].notna()].copy()

vv_species = sorted(df_vv["fair_specie"].unique().tolist())
print(f"  Mapped {len(vv_species)} VV species to FaIR v1.4.5 names.")

# Pivot each scenario: (fair_specie, year) -> value
smith_species_list = df_smith_full["variable"].tolist()

# Build one wide-format DataFrame per scenario, save to disk, and stash for the FaIR fill.
scenario_dfs = {}
for scen in scenarios:
    sub = df_vv[df_vv["scenario"] == scen]
    # pivot to species x year
    piv = sub.pivot_table(
        index="fair_specie", columns="year", values="value", aggfunc="first",
    )
    # Ensure all VV years 2023..2100 present in the columns
    for y in range(SMITH_LAST_YEAR + 1, VV_LAST_YEAR + 1):
        if y not in piv.columns:
            piv[y] = np.nan
    piv = piv[list(range(SMITH_LAST_YEAR + 1, VV_LAST_YEAR + 1))]

    # Build the wide CSV in canonical Smith format
    out_rows = []
    species_not_in_vv = []
    for v in smith_species_list:
        smith_row = df_smith_full[df_smith_full["variable"] == v].iloc[0]
        row = {
            "scenario": scen,
            "region":   smith_row["region"],
            "variable": v,
            "unit":     smith_row["unit"],
        }
        # Smith 1750.5..2022.5
        for c in keep_smith:
            row[c] = float(smith_row[c])
        # VV 2023.5..2100.5
        if v in piv.index:
            vv_vals = piv.loc[v].values.astype(float)
            scale = VV_TO_SMITH_SCALE.get(v, 1.0)
            vv_vals = vv_vals * scale
            for i, y in enumerate(range(SMITH_LAST_YEAR + 1, VV_LAST_YEAR + 1)):
                row[f"{y}.5"] = float(vv_vals[i])
        else:
            species_not_in_vv.append(v)
            # Hold constant at Smith 2022.5 value across the VV block
            last_val = float(smith_row[keep_smith[-1]])
            for y in range(SMITH_LAST_YEAR + 1, VV_LAST_YEAR + 1):
                row[f"{y}.5"] = last_val
        # Tail: hold at 2100.5 value
        v2100 = row[f"{VV_LAST_YEAR}.5"]
        for y in range(VV_LAST_YEAR + 1, TAIL_END_YEAR + 1):
            row[f"{y}.5"] = float(v2100)
        out_rows.append(row)

    df_out = pd.DataFrame(
        out_rows,
        columns=["scenario", "region", "variable", "unit", *all_year_cols],
    )
    scenario_dfs[scen] = df_out

    slug = slugify(scen)
    out_csv = SPLICED_DIR / f"{slug}.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"    Wrote {out_csv.name}  "
          f"(species held constant from Smith: {len(species_not_in_vv)})")
