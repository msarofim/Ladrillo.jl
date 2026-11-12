#!/usr/bin/env python3
"""
build_region_inventory.py -- the per-RGI-region glacier inventory that
ladrillo_data.py splits into reservoirs.

PORT of Part 2 of the research script `python/diag_constraint_anatomy.py`
(SLR-RFF-BRICK @ 3f7f2aa, lines 203-265). Part 1 of that script (the
one-reservoir arithmetic for the 2026-08-07 constraint memo) and its `exec` of
`d0_glacier_shootout.py` are not ported: nothing in Part 2 reads them, and the
only output Ladrillo needs is the regions table. The numerical code below is
copied verbatim; the only changes are I/O (paths from --sources / --work, no
eras/corners tables, no figure).

Reads (under --sources):
  gmip3/table_S1a.csv           GlacierMIP3 (Zekollari et al. 2025) Table S1a:
                                committed loss at +1.2 K, response time at ~1.5 K
  gmip3/table_S3.csv            Table S3: glacier mass in 2020 (Gt)
  gmip3/3_shift_summary_region_characteristicsFeb12_2024.csv
                                response time at +3.0 K, 2000-2019 regional
                                warming, Hugonnet 2000-2019 volume change
Writes (under --work):
  outputs/diag_constraint_anatomy_regions.csv   (the research file name, %.4f)

Usage: python build_region_inventory.py --sources DIR --work DIR
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

OUT_REL = "outputs/diag_constraint_anatomy_regions.csv"
FLOAT_FORMAT = "%.4f"          # the research file's format
EXCLUDED_REGION = "05"         # Greenland periphery: in the GIS target, not the glacier one


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sources", required=True, help="data/sources directory")
    ap.add_argument("--work", required=True, help="build directory (outputs/ is written under it)")
    a = ap.parse_args()
    GMIP3_DIR = os.path.join(a.sources, "gmip3")
    OUT_REGS = os.path.join(a.work, OUT_REL)
    os.makedirs(os.path.dirname(OUT_REGS), exist_ok=True)

    # ---------------------------------------------------------------- Part 2: regions
    s1a = pd.read_csv(os.path.join(GMIP3_DIR, "table_S1a.csv"), header=None, skiprows=3)
    s3 = pd.read_csv(os.path.join(GMIP3_DIR, "table_S3.csv"), index_col=0)
    s3 = s3[s3["rgi_reg"] != "global"]
    regchar = pd.read_csv(os.path.join(GMIP3_DIR,
                                       "3_shift_summary_region_characteristicsFeb12_2024.csv"),
                          index_col=0)
    regchar = regchar[regchar["rgi_reg"] != "All"]
    regchar["reg"] = regchar["rgi_reg"].map(lambda r: f"{int(r):02d}")
    regchar = regchar.set_index("reg")

    mass = {f"{int(r):02d}": m for r, m in
            zip(s3["rgi_reg"], s3["Glacier mass in 2020$^b$ (Gt)"])}

    def parse_committed(cell):
        m = re.match(r"\s*(-?\d+)\s*\[", str(cell))
        return float(m.group(1)) if m else np.nan

    def parse_resp(cell):
        m = re.match(r"\s*(-?\d+)\s*\[", str(cell))
        return float(m.group(1)) if m else np.nan

    reg_rows = []
    for _, row in s1a.iterrows():
        reg = str(row.iloc[10]).strip()
        if reg in ("global", "nan"):
            continue
        reg = f"{int(reg):02d}"
        if reg == EXCLUDED_REGION:
            continue                                # excluded from our scope
        reg_rows.append(dict(
            reg=reg,
            mass_gt=mass[reg],
            committed12_pct=parse_committed(row.iloc[2]),
            resp_time_15C_yr=parse_resp(row.iloc[9]),
            resp_time_30C_yr=float(regchar.loc[reg, "resp_time_-50%_3_0_deg"]),
            warming_2000_19_K=float(regchar.loc[reg, "temp_ch_avg_2000-2019_vs_1901-1920"]),
            melt_rate_gt_yr=-float(regchar.loc[reg, "dvoldt_m3_hugonnet"]) * 900.0 / 1e12,
        ))
    regs = pd.DataFrame(reg_rows)
    regs["stock_share"] = regs.mass_gt / regs.mass_gt.sum()
    regs["committed_share"] = (regs.committed12_pct / 100 * regs.mass_gt)
    regs["committed_share"] /= regs["committed_share"].sum()
    regs["melt_share"] = regs.melt_rate_gt_yr / regs.melt_rate_gt_yr.sum()
    regs = regs.sort_values("mass_gt", ascending=False)

    regs.to_csv(OUT_REGS, index=False, float_format=FLOAT_FORMAT)
    print(f"build_region_inventory | {len(regs)} regions (excl r{EXCLUDED_REGION}) | wrote {OUT_REL}")


if __name__ == "__main__":
    main()
