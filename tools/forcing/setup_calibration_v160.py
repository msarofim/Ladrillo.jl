"""
setup_calibration_v160.py

Assemble a `calibration_v160/` directory laid out exactly like `calibration_v145/`
so the existing FaIR drivers can run fair-calibrate **v1.6.0** by pointing
`--calib-dir` at it.

SOURCE: Zenodo record 18828694, "fair calibration data" v1.6.0, 2026-03-02
(fair-calibrate release v1.6.0, "FASTMIP Phase 1", constraint set `fastmip v1`).
Individual files are fetched; the 2.9 GB `fair_calibrate.zip` is not needed.

WHAT 1.6.0 CHANGES STRUCTURALLY vs 1.4.5 (this is why a plain parameter swap is
not enough):
  1. A new species **Irrigation** (62 species vs 61) with its own
     `forcing_scale[Irrigation]`, driven as PRESCRIBED forcing.
  2. **Land use** moves from `input_mode=calculated` to `forcing`, shipped as a
     per-ScenarioMIP-scenario timeseries.
  3. Volcanic and solar ship as separate CMIP7 files rather than one combined
     `volcanic_solar.csv`.

TWO DELIBERATE DEVIATIONS, made so that a van-Vuuren run and an RFF-SP run are
COMPARABLE (the whole point is a percentile placement of markers inside the RFF
band, which is meaningless if the two sides get different forcing treatments):

  * **Land use is set back to `calculated`.** 1.6.0 ships land-use forcing only
    for the seven ScenarioMIP markers; RFF-SP has no such file. `calculated`
    derives it from cumulative CO2 AFOLU via
    `land_use_cumulative_emissions_to_forcing` (1.6.0 value -0.00021836, vs
    1.4.5's -0.00022247, a 1.8% difference), which is what the v1.6.0 release
    note describes for the future anyway ("Land use forcing in the future
    follows cumulative land use CO2 emissions"). This keeps land use
    scenario-responsive on BOTH sides.

  * **Irrigation is applied as a single SHARED trajectory** (the mean across the
    seven markers) to every scenario on both sides. Irrigation is driven by an
    ISIMIP/population dataset, not by our emissions, and has no RFF-SP analogue.
    It is small and nearly common-mode (2100: -0.053 to -0.087 W/m2, spread
    0.034), so as a shared offset it cancels almost exactly in a percentile
    placement. Using per-marker irrigation for VV while giving RFF none would
    introduce exactly the asymmetry the test must avoid.

Both deviations are logged at build time and recorded in the output README.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/setup_calibration_v160.py @ 8984a13. The assembly is copied UNCHANGED;
only the I/O plumbing differs:
  * the five Zenodo files are READ from --raw-dir instead of fetched over the
    network. The download lives in fetch_inputs.py, which writes them to that
    directory and checks each against inputs.sha256. The bytes are parsed
    through io.BytesIO exactly as before, so a file on disk and the same file
    fetched give the same result.
  * --out-dir is REQUIRED (was FaIRtoFrEDI/calibration_v160).
Normally driven by build_forcing.py. By hand:
    python setup_calibration_v160.py --raw-dir INPUTS/zenodo_18828694 \
        --out-dir WORK/calibration_v160
"""

from pathlib import Path
import argparse
import io
import sys

import numpy as np
import pandas as pd

ZENODO_RECORD = "18828694"
ZENODO_LABEL  = ("Zenodo 10.5281/zenodo.18828694 'fair calibration data' v1.6.0 "
                 "(2026-03-02); fair-calibrate v1.6.0 'FASTMIP Phase 1', "
                 "constraint set 'fastmip v1'")

NEED = [
    "calibrated_constrained_parameters.csv",
    "species_configs_properties.csv",
    "volcanic_forcing_timebounds_cmip7.csv",
    "solar_forcing_timebounds_cmip7.csv",
    "irrigation_forcing_timebounds_cmip7.csv",
]

# Year grid must match the drivers' define_time(START_YEAR, END_YEAR) plus one,
# mirroring calibration_v145/volcanic_solar.csv which spans 1750..2302.
FORCING_FIRST_YEAR = 1750
FORCING_LAST_YEAR  = 2302

LANDUSE_MODE_FOR_SYMMETRY = "calculated"   # see module docstring


def read_raw(raw_dir, name):
    """The bytes fetch() used to return, read from a local copy instead."""
    p = Path(raw_dir) / name
    if not p.exists():
        sys.exit(f"ERROR: missing {p} -- run fetch_inputs.py first")
    print(f"  reading {name} ...")
    return p.read_bytes()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", type=Path, required=True,
                    help=f"directory holding the Zenodo {ZENODO_RECORD} files")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    OUT_DIR = args.out_dir
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("── assembling calibration_v160/ ──────────────────────────────────────")
    print(f"  source: {ZENODO_LABEL}")

    raw = {n: read_raw(args.raw_dir, n) for n in NEED}

    # ── parameters ───────────────────────────────────────────────────────────
    params = pd.read_csv(io.BytesIO(raw["calibrated_constrained_parameters.csv"]),
                         index_col=0)
    params.to_csv(OUT_DIR / "calibrated_constrained_parameters_1.6.0.csv")
    print(f"\n  parameters: {params.shape[0]} configs x {params.shape[1]} cols")

    # ── species configs, with Land use returned to `calculated` ─────────────
    spc = pd.read_csv(io.BytesIO(raw["species_configs_properties.csv"]))
    lu = spc["name"] == "Land use"
    before = spc.loc[lu, "input_mode"].iloc[0]
    spc.loc[lu, "input_mode"] = LANDUSE_MODE_FOR_SYMMETRY
    spc.to_csv(OUT_DIR / "species_configs_properties_1.6.0.csv", index=False)
    print(f"  species: {len(spc)}  |  Land use input_mode {before!r} -> "
          f"{LANDUSE_MODE_FOR_SYMMETRY!r} (DEVIATION, for VV/RFF symmetry)")
    print(f"  land_use_cumulative_emissions_to_forcing = "
          f"{float(spc.loc[spc.name=='CO2 AFOLU','land_use_cumulative_emissions_to_forcing'].iloc[0]):.8f}")

    # ── combined forcing CSV: Volcanic + Solar + (shared) Irrigation ────────
    volc = pd.read_csv(io.BytesIO(raw["volcanic_forcing_timebounds_cmip7.csv"]),
                       index_col=0).iloc[:, 0]
    solr = pd.read_csv(io.BytesIO(raw["solar_forcing_timebounds_cmip7.csv"]),
                       index_col=0).iloc[:, 0]
    irr = pd.read_csv(io.BytesIO(raw["irrigation_forcing_timebounds_cmip7.csv"]),
                      index_col=0)
    irr_shared = irr.mean(axis=1)
    print(f"\n  irrigation: shared = mean of {list(irr.columns)}")
    for y in (2100, 2300):
        if y in irr.index:
            print(f"    {y}: per-marker {irr.loc[y].min():+.4f}..{irr.loc[y].max():+.4f}"
                  f"  -> shared {irr_shared.loc[y]:+.4f} W/m2  (DEVIATION)")

    years = list(range(FORCING_FIRST_YEAR, FORCING_LAST_YEAR + 1))

    def on_grid(s, label):
        s = s[~s.index.duplicated()]
        out = s.reindex(years)
        if out.isna().any():
            # CMIP7 series end at 2501 so a gap can only be at the front;
            # hold the last available value forward, and report it.
            n_missing = int(out.isna().sum())
            out = out.ffill().bfill()
            print(f"    WARN {label}: {n_missing} year(s) filled by hold")
        return out.to_numpy(dtype=float)

    rows = [
        {"Scenario": "all", "Variable": "Volcanic", "Region": "World", "Unit": "W/m2",
         **dict(zip(map(str, years), on_grid(volc, "Volcanic")))},
        {"Scenario": "all", "Variable": "Solar", "Region": "World", "Unit": "W/m2",
         **dict(zip(map(str, years), on_grid(solr, "Solar")))},
        {"Scenario": "all", "Variable": "Irrigation", "Region": "World", "Unit": "W/m2",
         **dict(zip(map(str, years), on_grid(irr_shared, "Irrigation")))},
    ]
    fdf = pd.DataFrame(rows, columns=["Scenario", "Variable", "Region", "Unit",
                                      *map(str, years)])
    fdf.to_csv(OUT_DIR / "volcanic_solar.csv", index=False)
    print(f"  forcing CSV: {len(fdf)} rows x {len(years)} years "
          f"({years[0]}..{years[-1]})")

    (OUT_DIR / "README.md").write_text(f"""# calibration_v160

Assembled by `scripts/setup_calibration_v160.py`. Do not hand-edit.

**Source:** {ZENODO_LABEL}

Laid out like `calibration_v145/` so the existing drivers can use it via
`--calib-dir calibration_v160`.

## Files
- `calibrated_constrained_parameters_1.6.0.csv` — {params.shape[0]} configs,
  {params.shape[1]} columns (adds `forcing_scale[Irrigation]` vs 1.4.5).
- `species_configs_properties_1.6.0.csv` — {len(spc)} species (adds `Irrigation`).
- `volcanic_solar.csv` — Volcanic + Solar (CMIP7) + shared Irrigation,
  `Scenario="all"`, {years[0]}..{years[-1]}.

## Two deviations from stock 1.6.0 (both for VV-vs-RFF comparability)
1. **Land use** `input_mode` set {before!r} -> {LANDUSE_MODE_FOR_SYMMETRY!r}.
   1.6.0 ships land-use forcing only for the 7 ScenarioMIP markers; RFF-SP has
   no such file. `calculated` derives it from cumulative CO2 AFOLU, which is the
   rule 1.6.0 itself states for the future, and keeps land use
   scenario-responsive on both sides.
2. **Irrigation** applied as ONE shared trajectory (mean of the 7 markers) to
   every scenario. It is ISIMIP/population-driven with no RFF-SP analogue, is
   small and nearly common-mode, and as a shared offset cancels in a percentile
   placement.

Neither deviation is appropriate if you want to reproduce van Vuuren's own
figure exactly; both are appropriate for comparing markers against the RFF-SP
distribution on equal footing.
""")
    print(f"\nWrote {OUT_DIR}")


if __name__ == "__main__":
    main()
