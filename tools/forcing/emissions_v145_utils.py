"""
emissions_v145_utils.py -- the subset of the research repo's
FaIRtoFrEDI/emissions_v145_utils.py @ 8984a13 that build_vanvuuren_spliced.py
needs: the RFF-SP -> FaIR species map, the Mt -> Gt CO2 scale, the Smith
last-year constant, and load_smith(). Copied as source text (extracted with
ast.get_source_segment, not retyped).

PORT: load_smith() reads SMITH_FILE, which the research module hard-coded to
calibration_v145/smith2024_extended_1750-2025.csv. It is now set by the caller
(build_vanvuuren_spliced.py --smith-file) before load_smith() is called. That
file is chrisroadmap/temperature-attribution output/emissions_extrapolated_1750-2025.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd

SMITH_FILE = None      # set by the caller; was Path("calibration_v145/smith2024_extended_1750-2025.csv")

SMITH_LAST_YEAR_INT = 2020
PREFIX = "AR6 climate diagnostics|Emissions|"
SPECIES_MAP = {
    "BC": "BC", "OC": "OC", "NH3": "NH3", "NOx": "NOx", "VOC": "VOC",
    "CO": "CO", "Sulfur": "Sulfur", "CH4": "CH4", "N2O": "N2O",
    "NF3": "NF3", "SF6": "SF6", "SO2F2": "SO2F2",
    "CO2|Energy and Industrial Processes": "CO2 FFI",
    "CO2|AFOLU": "CO2 AFOLU",
    "CCl4": "CCl4", "CH2Cl2": "CH2Cl2", "CHCl3": "CHCl3",
    "CH3Br": "CH3Br", "CH3CCl3": "CH3CCl3", "CH3Cl": "CH3Cl",
    "CFC11": "CFC-11", "CFC12": "CFC-12", "CFC113": "CFC-113",
    "CFC114": "CFC-114", "CFC115": "CFC-115",
    "HCFC22": "HCFC-22", "HCFC141b": "HCFC-141b", "HCFC142b": "HCFC-142b",
    "Halon1211": "Halon-1211", "Halon1301": "Halon-1301",
    "Halon2402": "Halon-2402",
    "HFC|HFC23": "HFC-23", "HFC|HFC32": "HFC-32", "HFC|HFC125": "HFC-125",
    "HFC|HFC134a": "HFC-134a", "HFC|HFC143a": "HFC-143a",
    "HFC|HFC152a": "HFC-152a", "HFC|HFC227ea": "HFC-227ea",
    "HFC|HFC236fa": "HFC-236fa", "HFC|HFC245fa": "HFC-245fa",
    "HFC|HFC365mfc": "HFC-365mfc", "HFC|HFC43-10": "HFC-4310mee",
    "PFC|CF4": "CF4", "PFC|C2F6": "C2F6", "PFC|C3F8": "C3F8",
    "PFC|C4F10": "C4F10", "PFC|C5F12": "C5F12", "PFC|C6F14": "C6F14",
    "PFC|C7F16": "C7F16", "PFC|C8F18": "C8F18", "PFC|cC4F8": "c-C4F8",
}
RFF_TO_SMITH_SCALE = {"CO2 FFI": 1e-3, "CO2 AFOLU": 1e-3}


def load_smith():
    """Load Smith 2024 historical, returning DataFrame and the years to keep."""
    df_smith = pd.read_csv(SMITH_FILE)
    smith_year_cols = [c for c in df_smith.columns
                       if c not in ("scenario", "region", "variable", "unit")]
    keep_smith = [c for c in smith_year_cols
                  if float(c) <= SMITH_LAST_YEAR_INT + 0.5]
    return df_smith, keep_smith
