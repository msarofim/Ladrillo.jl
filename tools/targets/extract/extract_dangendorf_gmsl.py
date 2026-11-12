#!/usr/bin/env python3
"""
extract_dangendorf_gmsl.py -- MAINTAINER STEP: annual global-mean sea level from
the Dangendorf et al. (2024) Kalman-smoother fields, i.e.
data/sources/dangendorf2024/dangendorf2024_gmsl_annual.csv.

PORT of lines 83-94 of the research `python/diag_dangendorf_vs_frederikse.py`
(SLR-RFF-BRICK @ a954701), verbatim apart from I/O. The rest of that script
(the Dangendorf-vs-Frederikse comparison, its figure and summary) is not
needed and not ported.

Why the fields and not the published global file: the global file in Zenodo
10621070 is slot-shifted upstream (its GMSLHR holds the barystatic series).
The true GMSL is the cos(lat)-weighted mean of the "HR" field, as the record's
own Master_Final.m computes it. Dangendorf's corrected v2 global file
(2026-08-07, shipped in data/sources/dangendorf2024/) matches this extraction
to 0.0000 mm after removing a constant 53.45 mm baseline offset.
The `se_mm_indicative` column is read from the slot-shifted published file
and is NOT used by the targets (prep_recalib_targets_ext.py takes the SE from
the v2 file); it is kept so the file is byte-identical to the research one.

Units: m in the file, mm out. Years 1900-2021.

Usage: python extract_dangendorf_gmsl.py --raw-dir RAW --out-dir DIR
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import target_inputs as TI  # noqa: E402

FIELDS_REL = "dangendorf2024_KalmanSmootherHR_Fields.nc"
GLOBAL_REL = "dangendorf2024_KalmanSmootherHR_Global.nc"
GRID_REL = "dangendorf2024_grid_lonlat.csv"
OUT_NAME = "dangendorf2024_gmsl_annual.csv"
YEARS = (1900, 2021)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    pins = TI.read_pins()
    for rel in (FIELDS_REL, GLOBAL_REL, GRID_REL):
        got = TI.sha256(os.path.join(a.raw_dir, rel))
        if got != pins["raw/" + rel]:
            sys.exit(f"ABORT: {rel} sha256 {got} != pinned {pins['raw/' + rel]} (inputs.sha256)")

    # ---- verbatim (diag_dangendorf_vs_frederikse.py:87-94) ----
    grid = pd.read_csv(os.path.join(a.raw_dir, GRID_REL))
    wlat = np.cos(np.deg2rad(grid["lat"].values))
    fields = xr.open_dataset(os.path.join(a.raw_dir, FIELDS_REL), decode_times=False)
    d_years = np.arange(YEARS[0], YEARS[1] + 1)
    d_gmsl = (np.asarray(fields["HR"]) @ wlat) / wlat.sum() * 1000.0
    d_se = np.asarray(xr.open_dataset(os.path.join(a.raw_dir, GLOBAL_REL))["GMSLHRSE"])[0] * 1000.0
    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, OUT_NAME)
    pd.DataFrame({"year": d_years, "gmsl_mm": d_gmsl, "se_mm_indicative": d_se}).to_csv(out, index=False)
    print(f"extract_dangendorf_gmsl | {YEARS[0]}-{YEARS[1]} | wrote {out}")


if __name__ == "__main__":
    main()
