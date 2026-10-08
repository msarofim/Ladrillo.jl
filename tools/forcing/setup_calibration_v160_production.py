"""
setup_calibration_v160_production.py

Assemble `calibration_v160_prod/` -- the PRODUCTION fair-calibrate 1.6.0 stack,
as distinct from `calibration_v160/`, which exists for the van-Vuuren-vs-RFF
percentile comparison and carries two deliberate deviations for that purpose.

WHY A SECOND DIRECTORY RATHER THAN EDITING THE FIRST
`calibration_v160/` forces Land use back to `input_mode=calculated` and applies
Irrigation as one shared 7-marker mean. Both are CORRECT for a percentile
placement of markers inside an RFF-SP band, which is meaningless if the two
sides get different forcing treatments, and both are WRONG for production:
`diag_calib160_igcc_ladder.py` measured the land-use deviation at 0.025 K of
historical GMST fit (rung R5 -> R6). Marcus 2026-08-28 chose the hybrid --
prescribed 1.6.0 forcing wherever a scenario maps to a CMIP7 marker, the
`calculated`/shared treatment only where none does (RFF-SP). Keeping both
directories means neither use silently inherits the other's choice.

WHAT THIS DIRECTORY CONTAINS
  - the same 841-config parameter file as `calibration_v160/`
  - species configs with Land use at STOCK 1.6.0 `input_mode=forcing`
  - one prescribed-forcing CSV PER CMIP7 MARKER, each carrying Volcanic, Solar,
    Land use and Irrigation, so a scenario picks its marker and gets all four

THE MARKER MAPPING IS THE OPEN PIECE. 1.6.0 ships land-use and irrigation
forcing for the seven CMIP7 ScenarioMIP markers; our drivers run CMIP6 SSPs, and
no CMIP7 marker IS a CMIP6 SSP. Measured here at build time: across the seven
markers the two series are BIT-IDENTICAL through 2023 and diverge only after,
reaching 0.0895 (land use) and 0.0342 (irrigation) W/m2 by 2100. So the mapping
cannot affect any historical result, and this script writes every marker rather
than choosing for the caller. Drivers name their marker explicitly.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/setup_calibration_v160_production.py @ 8984a13. The assembly is copied
UNCHANGED; only the I/O plumbing differs: the three directories that were
FaIRtoFrEDI/{calibration_v160, v160_cache, calibration_v160_prod} are the
REQUIRED flags --src-calib-dir, --raw-dir (the Zenodo 18828694 files) and
--out-dir.
Normally driven by build_forcing.py. By hand:
    python setup_calibration_v160_production.py \
        --src-calib-dir WORK/calibration_v160 --raw-dir INPUTS/zenodo_18828694 \
        --out-dir WORK/calibration_v160_prod
"""

from pathlib import Path
import argparse
import shutil
import sys

import numpy as np
import pandas as pd

# ── Named constants ──────────────────────────────────────────────────────────
CALIB_VERSION = "1.6.0"
SOURCE_LABEL = ("Zenodo 10.5281/zenodo.18828694 'fair calibration data' v1.6.0 "
                "(2026-03-02); fair-calibrate v1.6.0 'FASTMIP Phase 1', "
                "constraint set 'fastmip v1'")

# Filenames; their directories come from the CLI (see main()).
PARAMS_NAME = f"calibrated_constrained_parameters_{CALIB_VERSION}.csv"
SPECIES_NAME = f"species_configs_properties_{CALIB_VERSION}.csv"
VOLCANIC_NAME = "volcanic_forcing_timebounds_cmip7.csv"
SOLAR_NAME = "solar_forcing_timebounds_cmip7.csv"
LANDUSE_NAME = "land_use_forcing_timebounds_cmip7.csv"
IRRIGATION_NAME = "irrigation_forcing_timebounds_cmip7.csv"

# Year grid mirrors calibration_v145/volcanic_solar.csv (1750..2302), one past
# the drivers' END_YEAR so define_time never runs off the end of the forcing.
FORCING_FIRST_YEAR = 1750
FORCING_LAST_YEAR = 2302

# The two scenario-resolved species. Named here so the label in the README, the
# divergence report and the per-marker CSVs all derive from one list. Filled
# with full paths in main(), once --raw-dir is known.
PER_MARKER_SPECIES = {"Land use": LANDUSE_NAME, "Irrigation": IRRIGATION_NAME}

# Land use is PRESCRIBED in stock 1.6.0; calibration_v160/ overrode it. Assert we
# are undoing exactly that override rather than assuming.
LANDUSE_MODE_PRODUCTION = "forcing"


def on_grid(s, years, label):
    s = s[~s.index.duplicated()]
    out = s.reindex(years)
    if out.isna().any():
        n = int(out.isna().sum())
        out = out.ffill().bfill()
        print(f"    WARN {label}: {n} year(s) filled by hold")
    return out.to_numpy(float)


def main():
    global SRC_CALIB, OUT_DIR, PARAMS_SRC, SPECIES_SRC, VOLCANIC_SRC, SOLAR_SRC
    global LANDUSE_SRC, IRRIGATION_SRC, PER_MARKER_SPECIES
    ap = argparse.ArgumentParser()
    ap.add_argument("--src-calib-dir", type=Path, required=True,
                    help="calibration_v160/ (params + species come from here)")
    ap.add_argument("--raw-dir", type=Path, required=True,
                    help="the Zenodo 18828694 files (raw CMIP7 forcing)")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    SRC_CALIB, V160_CACHE, OUT_DIR = args.src_calib_dir, args.raw_dir, args.out_dir
    PARAMS_SRC = SRC_CALIB / PARAMS_NAME
    SPECIES_SRC = SRC_CALIB / SPECIES_NAME
    VOLCANIC_SRC = V160_CACHE / VOLCANIC_NAME
    SOLAR_SRC = V160_CACHE / SOLAR_NAME
    LANDUSE_SRC = V160_CACHE / LANDUSE_NAME
    IRRIGATION_SRC = V160_CACHE / IRRIGATION_NAME
    PER_MARKER_SPECIES = {k: V160_CACHE / v for k, v in PER_MARKER_SPECIES.items()}

    for p in (PARAMS_SRC, SPECIES_SRC, VOLCANIC_SRC, SOLAR_SRC,
              LANDUSE_SRC, IRRIGATION_SRC):
        if not p.exists():
            sys.exit(f"ERROR: missing input {p}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("── assembling calibration_v160_prod/ ────────────────────────────────")
    print(f"  source: {SOURCE_LABEL}")

    # ── parameters: byte-identical copy ──────────────────────────────────────
    dst_params = OUT_DIR / PARAMS_SRC.name
    shutil.copyfile(PARAMS_SRC, dst_params)
    params = pd.read_csv(dst_params, index_col=0)
    print(f"\n  parameters: {params.shape[0]} configs x {params.shape[1]} cols "
          f"(copied verbatim from {SRC_CALIB.name}/)")

    # ── species configs: Land use back to STOCK 1.6.0 ────────────────────────
    spc = pd.read_csv(SPECIES_SRC)
    lu = spc["name"] == "Land use"
    before = spc.loc[lu, "input_mode"].iloc[0]
    if before == LANDUSE_MODE_PRODUCTION:
        sys.exit(f"ERROR: {SPECIES_SRC.name} already has Land use input_mode "
                 f"{before!r}; this script expects to be undoing the "
                 f"calibration_v160 'calculated' deviation. Check the source.")
    spc.loc[lu, "input_mode"] = LANDUSE_MODE_PRODUCTION
    dst_species = OUT_DIR / SPECIES_SRC.name
    spc.to_csv(dst_species, index=False)
    print(f"  species: {len(spc)}  |  Land use input_mode {before!r} -> "
          f"{LANDUSE_MODE_PRODUCTION!r} (= stock 1.6.0; reverts the "
          f"{SRC_CALIB.name} deviation)")

    # ── per-marker prescribed forcing ────────────────────────────────────────
    years = list(range(FORCING_FIRST_YEAR, FORCING_LAST_YEAR + 1))
    volc = pd.read_csv(VOLCANIC_SRC, index_col=0).iloc[:, 0]
    solr = pd.read_csv(SOLAR_SRC, index_col=0).iloc[:, 0]
    per_marker = {k: pd.read_csv(v, index_col=0) for k, v in PER_MARKER_SPECIES.items()}

    markers = list(per_marker["Land use"].columns)
    for name, df in per_marker.items():
        if list(df.columns) != markers:
            sys.exit(f"ERROR: {name} markers {list(df.columns)} != {markers}")
    print(f"\n  markers: {markers}")

    # Report WHERE the marker choice can and cannot matter, measured not assumed.
    print("  across-marker spread (W/m2) -- the mapping question's whole size:")
    for name, df in per_marker.items():
        row = f"    {name:<11}"
        for y in (1900, 2000, 2023, 2050, 2100, 2300):
            if y in df.index:
                row += f"  {y}: {float(df.loc[y].max() - df.loc[y].min()):.4f}"
        print(row)
    hist = per_marker["Land use"].loc[:2023]
    identical = bool(np.allclose(hist.to_numpy(), hist.iloc[:, [0]].to_numpy(), atol=0, rtol=0))
    print(f"    land use bit-identical across markers through 2023: {identical}")
    if not identical:
        print("    NOTE: NOT bit-identical -- the mapping now touches history too, "
              "which invalidates the 'history is marker-free' argument.")

    for m in markers:
        rows = [
            {"Scenario": "all", "Variable": "Volcanic", "Region": "World", "Unit": "W/m2",
             **dict(zip(map(str, years), on_grid(volc, years, "Volcanic")))},
            {"Scenario": "all", "Variable": "Solar", "Region": "World", "Unit": "W/m2",
             **dict(zip(map(str, years), on_grid(solr, years, "Solar")))},
        ]
        for name, df in per_marker.items():
            rows.append({"Scenario": "all", "Variable": name, "Region": "World",
                         "Unit": "W/m2",
                         **dict(zip(map(str, years), on_grid(df[m], years, f"{name}/{m}")))})
        out = pd.DataFrame(rows, columns=["Scenario", "Variable", "Region", "Unit",
                                          *map(str, years)])
        out.to_csv(OUT_DIR / f"volcanic_solar_{m}.csv", index=False)
    print(f"\n  wrote {len(markers)} per-marker forcing CSVs "
          f"(Volcanic + Solar + {' + '.join(PER_MARKER_SPECIES)}), "
          f"{years[0]}..{years[-1]}")

    (OUT_DIR / "README.md").write_text(f"""# calibration_v160_prod

Assembled by `tools/forcing/setup_calibration_v160_production.py`. Do not hand-edit.

**Source:** {SOURCE_LABEL}

The **production** fair-calibrate {CALIB_VERSION} stack. Distinct from
`calibration_v160/`, which exists for the van-Vuuren-vs-RFF percentile
comparison and deliberately deviates from stock 1.6.0 in two ways (Land use
forced to `calculated`; Irrigation as one shared 7-marker mean). Those
deviations are right for that comparison and wrong for production --
`diag_calib160_igcc_ladder.py` measured the land-use one at **0.025 K** of
historical GMST fit. Marcus chose the hybrid on 2026-08-28: prescribed 1.6.0
forcing wherever a scenario maps to a CMIP7 marker; `calculated` land use and
shared irrigation only where none does (RFF-SP).

## Files
- `{PARAMS_SRC.name}` — {params.shape[0]} configs x {params.shape[1]} columns,
  byte-identical to `calibration_v160/`.
- `{SPECIES_SRC.name}` — {len(spc)} species, Land use at **stock 1.6.0
  `input_mode=forcing`**.
- `volcanic_solar_<MARKER>.csv` for MARKER in {markers} — Volcanic + Solar
  (CMIP7) + per-marker Land use + per-marker Irrigation.

## The marker mapping
1.6.0 ships land-use and irrigation forcing for the seven CMIP7 ScenarioMIP
markers. Our drivers run CMIP6 SSPs, and no CMIP7 marker **is** a CMIP6 SSP, so
a driver must name the marker it uses. Measured at build time, across the seven
markers the two series are **bit-identical through 2023** and diverge only
after, reaching {float(per_marker['Land use'].loc[2100].max() - per_marker['Land use'].loc[2100].min()):.4f}
(land use) and {float(per_marker['Irrigation'].loc[2100].max() - per_marker['Irrigation'].loc[2100].min()):.4f}
(irrigation) W/m2 by 2100. **The mapping therefore cannot affect any historical
result** — including the Ladrillo calibration objective, whose 1850-2026 fit
window has only 3 years past the 2023 emissions splice, bounded there at about
0.002 W/m2. It does affect projections to 2300 and is an open question for
`ssp126`/`ssp585`.

This directory writes **every** marker and chooses for nobody.
""")
    print(f"\nWrote {OUT_DIR}")


if __name__ == "__main__":
    main()
