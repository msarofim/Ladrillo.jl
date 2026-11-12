#!/usr/bin/env python3
"""
build_lws_grace_extension.py -- extend the land-water (LWS) target past 2018
with GRACE/GRACE-FO.

PORT of the time-series half of the research
`python/build_lws_grace_extension.py --tag=L14` (SLR-RFF-BRICK @ d2e23d3). The
gridded half (masks, gates, the domain sum) is the maintainer step
extract/extract_grace_lws_domain.py, whose output is shipped as
data/sources/grace/grace_lws_domain_monthly.csv. From that series on, the code
is verbatim: deseasonalise (2 harmonics), annual means (>= 4 samples),
subtract GlaMBIE glaciers (17 regions: all but 0 = global, 5 = inside the
Greenland mask, 19 = no land cells at 0.5 degree), convert to cm SLE,
offset-match onto Frederikse's TWS over 2003-2018, hold 2024 at 2023.

THE CIRCULAR READ, REMOVED. The research script read the Frederikse TWS for
its splice offset from outputs/recalib_targets_ext.csv, the file that
prep_recalib_targets_ext.py writes AFTER reading this script's output. The two
only agreed because the 2003-2018 lws values in the targets file are
Frederikse's, whatever the extension. Here the Frederikse TWS is computed from
the Frederikse workbook with the same lines prep_recalib_targets_ext.py uses
(column mean / 10, minus its 1995-2005 mean), so the build order is acyclic:
this script, then prep. One detail of the old path is kept on purpose: the
values pass through an in-memory CSV write/read with pandas' defaults, because
that is how the research script received them and the default parser shifts
some of them by one ulp (see as_read_back_from_targets_csv). Without it the
extension moves by up to 1.1e-16 cm and the targets are no longer
byte-identical.

Dropped (nothing downstream reads them): the [RESULT] comparison with the old
hold-flat target, the console banner, the figure. `--tag` survives only as the
file-name convention the targets builder reads (L14; the extension does not
depend on a posterior).

Reads (under --sources): grace/grace_lws_domain_monthly.csv,
       glambie/calendar_years/*.csv, frederikse2020/frederikse2020_global_basin_timeseries.xlsx
Writes (under --work):   outputs/lws_grace_extension_L14.csv

Usage: python build_lws_grace_extension.py --sources DIR --work DIR
"""
import argparse
import glob
import io
import os

import numpy as np
import pandas as pd

TAG = "L14"
DOMAIN_REL = "grace/grace_lws_domain_monthly.csv"
GLAMBIE_REL = "glambie/calendar_years"
FRED_REL = "frederikse2020/frederikse2020_global_basin_timeseries.xlsx"
FRED_TWS_COL = "Terrestrial Water Storage [mean]"
OUT_REL = f"outputs/lws_grace_extension_{TAG}.csv"

GT_PER_CM_SLE = 3620.0     # 362 Gt/mm, the repo convention (prep_recalib_targets_ext.py)
GLACIER_EXCLUDE = {"5", "19"}   # 5 is inside the GIS mask; 19 has no land cells at 0.5 deg
N_HARMONICS = 2            # annual + semi-annual, for [SEAS]
MIN_MONTHS = 4             # a year needs this many samples even after deseasonalising
OVERLAP = (2003, 2018)     # offset-match window; Frederikse's TWS ends 2018
SPLICE_FROM = 2019
GLAMBIE_LAST = 2023        # GlaMBIE ends here -> see HOLD_LAST
HOLD_LAST = 2024           # Marcus 2026-08-24: hold 2024 flat from 2023 rather than drop it
EXT_LAST = 2026            # the target file's own horizon
BASE_Y0, BASE_Y1 = 1995, 2005   # the targets' reference window (prep_recalib_targets_ext.py)


def frederikse_tws(sources):
    """Frederikse 2020 TWS, cm, re-referenced to 1995-2005 with the mean's single offset:
    the same lines prep_recalib_targets_ext.py uses for every component (its `fred[tgt]`)."""
    g = pd.read_excel(os.path.join(sources, FRED_REL), "Global").rename(
        columns={"Unnamed: 0": "year"}).set_index("year")
    off = (g[FRED_TWS_COL] / 10.0).loc[BASE_Y0:BASE_Y1].mean()
    return (g[FRED_TWS_COL] / 10.0) - off


def as_read_back_from_targets_csv(s):
    """The values `s` as the research script SAW them: prep_recalib_targets_ext.py wrote them to
    recalib_targets_ext.csv with pandas' default float output (exact, shortest round-trip), and this
    script read them back with pandas' default parser, which mis-rounds some 17-digit decimals by one
    ulp (measured 2026-11-12, pandas 2.3.3: 5.6e-17 on the 2003-2018 TWS; float_precision="round_trip"
    gives 0). The splice offset, and so every extension value, carries that ulp. Sending the series
    through the same write/read pair, in memory, reproduces it without reading the targets file."""
    buf = io.StringIO()
    pd.DataFrame({"year": s.index, "lws": s.to_numpy()}).to_csv(buf, index=False)
    buf.seek(0)
    return pd.read_csv(buf).set_index("year")["lws"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sources", required=True)
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    OUT_CSV = os.path.join(a.work, OUT_REL)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

    dm = pd.read_csv(os.path.join(a.sources, DOMAIN_REL), float_precision="round_trip")
    t = pd.DatetimeIndex(pd.to_datetime(dm["time"]))
    g = dm["domain_mass_gt"].to_numpy(dtype=float)

    # ---- verbatim from here (build_lws_grace_extension.py:193-244) ----
    frac = (t.dayofyear - 1) / 365.25
    X = np.column_stack([np.ones(len(t))] +
                        [f(2 * np.pi * k * frac) for k in range(1, N_HARMONICS + 1)
                         for f in (np.sin, np.cos)])
    seas = X[:, 1:] @ np.linalg.lstsq(X, g, rcond=None)[0][1:]
    gd = pd.Series(g - seas, index=t.year)
    print(f"[SEAS]  removed a {seas.max() - seas.min():.0f} Gt peak-to-peak annual cycle "
          f"({N_HARMONICS} harmonics)")
    ann, nmo = gd.groupby(level=0).mean(), gd.groupby(level=0).count()
    yrs = np.array([y for y in range(2003, HOLD_LAST + 1)
                    if nmo.get(y, 0) >= MIN_MONTHS and y <= GLAMBIE_LAST])

    gl = {}
    for f in sorted(glob.glob(os.path.join(a.sources, GLAMBIE_REL, "*.csv"))):
        code = os.path.basename(f)[:-4].split("_")[0]
        if code == "0" or code in GLACIER_EXCLUDE:
            continue
        d = pd.read_csv(f)
        gl[code] = d.set_index(d.start_dates.astype(int))["combined_gt"]
    glac_cum = pd.DataFrame(gl).sum(axis=1).cumsum()
    print(f"[GLAC]  subtracting {len(gl)} GlaMBIE regions (all but "
          f"{sorted(GLACIER_EXCLUDE)}); {SPLICE_FROM}-{GLAMBIE_LAST} loss "
          f"{pd.DataFrame(gl).loc[SPLICE_FROM:GLAMBIE_LAST].sum().sum():.0f} Gt")

    lws = pd.Series(-(ann.loc[yrs].values - glac_cum.reindex(yrs).values) / GT_PER_CM_SLE,
                    index=yrs)

    # ---- splice onto the Frederikse frame (PORT: Frederikse TWS from the workbook, see the docstring)
    fred_lws = as_read_back_from_targets_csv(frederikse_tws(a.sources))
    ov = [y for y in yrs if OVERLAP[0] <= y <= OVERLAP[1]]
    off = fred_lws.loc[ov].mean() - lws.loc[ov].mean()
    lws = lws + off
    res = lws.loc[ov] - fred_lws.loc[ov]
    print(f"[OVERLAP] {ov[0]}-{ov[-1]} (n={len(ov)}), offset {off:+.3f} cm, rms "
          f"{np.sqrt((res ** 2).mean()):.3f}, max|res| {abs(res).max():.3f} cm")

    out = pd.Series(index=range(2003, EXT_LAST + 1), dtype=float)
    out.loc[lws.index] = lws.values
    for y in range(GLAMBIE_LAST + 1, HOLD_LAST + 1):
        out.loc[y] = out.loc[GLAMBIE_LAST]
    print(f"[HOLD]  {GLAMBIE_LAST + 1}-{HOLD_LAST} held flat at the {GLAMBIE_LAST} value "
          f"({out.loc[GLAMBIE_LAST]:.3f} cm)")

    pd.DataFrame({"year": out.index, "lws_grace_cm": out.values,
                  "is_held_flat": [(y > GLAMBIE_LAST) for y in out.index]}).to_csv(OUT_CSV, index=False)
    print(f"build_lws_grace_extension | wrote {OUT_REL}")


if __name__ == "__main__":
    main()
