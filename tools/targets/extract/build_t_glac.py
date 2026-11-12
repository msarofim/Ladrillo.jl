#!/usr/bin/env python3
"""
build_t_glac.py — glacier-area-weighted observed temperature series T_glac(t).

Option D data-prep (handoff_2026-08-06_glacier_structural_nu_decision.md §2):
the decadal glacier-melt shape is not a GMST shape (ETCW: Arctic warming
3-4x global by 1930-40 on most of the world's glacier area — memo §3d), so
we build the temperature the GLACIERS saw, to drive the Mengel/Nauels glacier
component historically in the D0 shootout.

Construction (all choices explicit; see PROVENANCE block written alongside):
  * Temperature field: HadCRUT5.0.2.0 ANALYSIS (statistically infilled)
    ensemble-mean monthly anomalies, 5-deg grid, 1850-present, rel 1961-1990.
    (Analysis variant chosen for pre-1930 Arctic coverage; the infilling is
    heaviest exactly in the early-Arctic window — documented caveat.)
  * Regions: GTN-G Glacier Regions 2023 first-order polygons
    (DOI 10.5904/gtng-glacreg-2023-07). The 2023 scheme splits RGI6 region 19
    into 19 (Subantarctic islands) + 20 (Antarctic mainland periphery); both
    are merged back into "r19" to match the RGI6 / GlaMBIE scope.
  * Region -> grid: fractional overlap of each 5-deg cell with the region
    polygon, via a SUBGRID x SUBGRID point-in-polygon test per cell
    (small regions like Iceland would otherwise catch zero cell centres).
  * Scope: regions 1-18 minus 5, plus 19 — the excl-r5-incl-r19 scope of the
    A2 inventory V = 0.290 (memo §2c; r5 lives in the GIS target).
  * Weights (headline): GlaMBIE per-region glacier area at year 2000
    (on-disk glambie_data.zip, DOI 10.5904/wgms-glambie-2024-07).
    Sensitivity: |cumulative 2000-2023 mass loss| share (melt-flux weights).
  * Annual = calendar-year mean of the 12 monthly values (per the
    Berkeley-parser calendar-year discipline); year kept only if >= MIN_MONTHS.
  * Baseline: anomalies re-referenced to 1850-1900 (multi-year baseline rule).

Outputs (under --out-dir; Ladrillo ships both as data/sources/hadcrut5/):
  t_glac_hadcrut5.csv          (year, t_glac_C, t_glac_massloss_C, gmst_hadcrut5_C, fill_frac)
  t_glac_regions_hadcrut5.csv  (year, r01..r19 annual anomalies rel 1850-1900)

MAINTAINER STEP. PORT of the research `python/build_t_glac.py` (SLR-RFF-BRICK
@ 597b7e1). The numerical code is verbatim. What changed (I/O only): the inputs
come from --raw-dir (HadCRUT5, the GTN-G shapefile) and --sources (the GlaMBIE
calendar-year files, which the research script read out of glambie_data.zip;
the bytes are the same members) and are checked against inputs.sha256 first;
the outputs go to --out-dir; the provenance markdown, the console report and
the figure are dropped.

Usage: python build_t_glac.py --raw-dir RAW --sources data/sources --out-dir DIR
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import target_inputs as TI  # noqa: E402

_ap = argparse.ArgumentParser(description="glacier-area-weighted HadCRUT5 per RGI region (maintainer step)")
_ap.add_argument("--raw-dir", required=True)
_ap.add_argument("--sources", required=True, help="data/sources (GlaMBIE calendar-year files)")
_ap.add_argument("--out-dir", required=True)
_ARGS = _ap.parse_args()
HAD_REL = "HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc"
SHP_REL = "GlacReg_2023/GTN-G_202307_o1regions.shp"
HAD_NC = os.path.join(_ARGS.raw_dir, HAD_REL)
SHP = os.path.join(_ARGS.raw_dir, SHP_REL)
GLAMBIE_DIR = os.path.join(_ARGS.sources, "glambie", "calendar_years")
OUT_CSV = os.path.join(_ARGS.out_dir, "t_glac_hadcrut5.csv")
OUT_REG = os.path.join(_ARGS.out_dir, "t_glac_regions_hadcrut5.csv")
_pins = TI.read_pins()
for _rel in (HAD_REL, SHP_REL, SHP_REL[:-4] + ".shx", SHP_REL[:-4] + ".dbf"):
    _got = TI.sha256(os.path.join(_ARGS.raw_dir, _rel))
    if _got != _pins["raw/" + _rel]:
        sys.exit(f"ABORT: {_rel} sha256 {_got} != pinned {_pins['raw/' + _rel]} (inputs.sha256)")

SCOPE_REGIONS = [r for r in range(1, 20) if r != 5]   # excl r5, incl r19 (V=0.290 scope)
MERGE_INTO_19 = {20}          # GTN-G 2023 splits RGI6 r19 into 19+20; merge back
SUBGRID = 5                   # per-5deg-cell fractional-overlap subsampling (SUBGRID^2 pts)
BASE0, BASE1 = 1850, 1900     # anomaly baseline
MIN_MONTHS = 12               # calendar-year completeness required
GLAMBIE_AREA_YEAR = 2000.0    # weight = glacier_area at this start_date

# ------------------------------------------------------------------ regions
import shapefile as pyshp                     # pyshp
from matplotlib.path import Path

sf = pyshp.Reader(SHP)
fields = [f[0] for f in sf.fields[1:]]
REG_KEY = "o1region" if "o1region" in fields else fields[0]

rings_by_region = {}
for sr in sf.shapeRecords():
    reg = int(sr.record[fields.index(REG_KEY)])
    if reg in MERGE_INTO_19:
        reg = 19
    shp = sr.shape
    parts = list(shp.parts) + [len(shp.points)]
    for i in range(len(shp.parts)):
        ring = np.asarray(shp.points[parts[i]:parts[i + 1]])
        if len(ring) >= 3:
            rings_by_region.setdefault(reg, []).append(Path(ring))

# ------------------------------------------------------------------ grid masks
import netCDF4

nc = netCDF4.Dataset(HAD_NC)
VAR = "tas_mean"
lat = nc.variables["latitude"][:]
lon = nc.variables["longitude"][:]
tv = nc.variables["time"]
dates = netCDF4.num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"))
yr_m = np.array([d.year for d in dates])
mo_m = np.array([d.month for d in dates])

dlat = np.abs(np.diff(lat)).mean()
dlon = np.abs(np.diff(lon)).mean()
off = (np.arange(SUBGRID) + 0.5) / SUBGRID - 0.5

frac = {}     # region -> (nlat, nlon) fractional overlap
for reg, rings in sorted(rings_by_region.items()):
    f = np.zeros((len(lat), len(lon)))
    for i, la in enumerate(lat):
        sub_la = la + off * dlat
        for j, lo in enumerate(lon):
            sub_lo = lo + off * dlon
            pts = np.array([(x, y) for y in sub_la for x in sub_lo])
            inside = np.zeros(len(pts), dtype=bool)
            for ring in rings:
                inside |= ring.contains_points(pts)
            # dateline wrap: retry shifted copies for polygons stored in 0..360 or crossing
            if not inside.any():
                for shift in (360.0, -360.0):
                    pts2 = pts.copy(); pts2[:, 0] += shift
                    for ring in rings:
                        inside |= ring.contains_points(pts2)
            f[i, j] = inside.mean()
    frac[reg] = f

# ------------------------------------------------------------------ weights
# PORT: the research read these members out of glambie_data.zip; same bytes, from data/sources/
names = {int(os.path.basename(n).split("_")[0]): n
         for n in sorted(glob.glob(os.path.join(GLAMBIE_DIR, "*.csv")))
         if os.path.basename(n).split("_")[0].isdigit()
         and int(os.path.basename(n).split("_")[0]) >= 1}
area_w, loss_w = {}, {}
for reg, n in names.items():
    df = pd.read_csv(n)
    row0 = df.loc[np.isclose(df["start_dates"], GLAMBIE_AREA_YEAR)]
    area_w[reg] = float(row0["glacier_area"].iloc[0])
    loss_w[reg] = abs(float(df["combined_gt"].sum()))

def normalize(w):
    tot = sum(w[r] for r in SCOPE_REGIONS)
    return {r: w[r] / tot for r in SCOPE_REGIONS}

W_AREA = normalize(area_w)
W_LOSS = normalize(loss_w)

# ------------------------------------------------------------------ per-region series
tas = nc.variables[VAR][:]                     # (time, lat, lon), masked array
tas = np.ma.filled(tas, np.nan)
coslat = np.cos(np.deg2rad(lat))[:, None]

reg_monthly = {}
cov_report = {}
for reg in sorted(frac):
    w2d = frac[reg] * coslat
    wsum_full = w2d.sum()
    if wsum_full == 0:
        raise RuntimeError(f"region {reg:02d}: polygon captured no grid area")
    wt = np.where(np.isfinite(tas), w2d[None, :, :], 0.0)
    num = np.nansum(np.where(np.isfinite(tas), tas, 0.0) * w2d[None, :, :], axis=(1, 2))
    den = wt.sum(axis=(1, 2))
    series = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
    reg_monthly[reg] = series
    cov_report[reg] = (float((frac[reg] > 0).sum()), float(den.min() / wsum_full),
                       float(den[yr_m < 1900].mean() / wsum_full) if (yr_m < 1900).any() else np.nan)

years = np.arange(yr_m.min(), yr_m.max() + 1)
def annualize(mseries):
    out = np.full(len(years), np.nan)
    for k, y in enumerate(years):
        sel = yr_m == y
        if sel.sum() >= MIN_MONTHS and np.isfinite(mseries[sel]).all():
            out[k] = mseries[sel].mean()
    return out

reg_annual = {reg: annualize(s) for reg, s in reg_monthly.items()}

# global HadCRUT5 (cos-lat mean of the full field) for amp_g reference
gnum = np.nansum(np.where(np.isfinite(tas), tas, 0.0) * coslat[None, :, :], axis=(1, 2))
gden = np.where(np.isfinite(tas), coslat[None, :, :], 0.0).sum(axis=(1, 2))
gmst_annual = annualize(gnum / np.maximum(gden, 1e-12))

# ---- fill missing region-YEARS (all pre-1898: early Arctic/HMA months with zero
# coverage) by per-region OLS on the global series in the NATIVE (1961-1990) frame:
# region ~ amp_r * global + c_r over all jointly-available years. Keeps every
# region complete back to 1850 so the composite and the 1850-1900 baseline exist;
# filled years carry scaled-global (no regional internal variability) — reported
# and stored in the fill_frac column.
fill_note = {}
fill_mask = {}
for reg, x in reg_annual.items():
    m = np.isfinite(x) & np.isfinite(gmst_annual)
    miss = ~np.isfinite(x) & np.isfinite(gmst_annual)
    fill_mask[reg] = miss
    if miss.any():
        A = np.vstack([gmst_annual[m], np.ones(m.sum())]).T
        amp_r, c_r = np.linalg.lstsq(A, x[m], rcond=None)[0]
        x[miss] = amp_r * gmst_annual[miss] + c_r
        fill_note[reg] = (int(miss.sum()), int(years[miss].min()), int(years[miss].max()),
                          float(amp_r))

def rebase(x):
    sel = (years >= BASE0) & (years <= BASE1)
    return x - np.nanmean(x[sel])

reg_annual = {r: rebase(x) for r, x in reg_annual.items()}
gmst_annual = rebase(gmst_annual)

t_glac_area = np.nansum([W_AREA[r] * reg_annual[r] for r in SCOPE_REGIONS], axis=0)
t_glac_loss = np.nansum([W_LOSS[r] * reg_annual[r] for r in SCOPE_REGIONS], axis=0)
ok = np.all([np.isfinite(reg_annual[r]) for r in SCOPE_REGIONS], axis=0) & np.isfinite(gmst_annual)
t_glac_area[~ok] = np.nan
t_glac_loss[~ok] = np.nan
fill_frac = np.sum([W_AREA[r] * fill_mask[r].astype(float) for r in SCOPE_REGIONS], axis=0)

# ------------------------------------------------------------------ outputs
os.makedirs(_ARGS.out_dir, exist_ok=True)
out = pd.DataFrame({"year": years[ok],
                    "t_glac_C": t_glac_area[ok],
                    "t_glac_massloss_C": t_glac_loss[ok],
                    "gmst_hadcrut5_C": gmst_annual[ok],
                    "fill_frac": fill_frac[ok]})
out.to_csv(OUT_CSV, index=False, float_format="%.5f")

reg_out = pd.DataFrame({"year": years})
for reg in sorted(reg_annual):
    reg_out[f"r{reg:02d}"] = reg_annual[reg]
reg_out = reg_out[np.isfinite(reg_out.drop(columns="year")).any(axis=1).to_numpy()]
reg_out.to_csv(OUT_REG, index=False, float_format="%.5f")

print(f"build_t_glac | wrote {OUT_CSV}, {OUT_REG}")
