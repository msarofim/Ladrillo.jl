#!/usr/bin/env python3
"""
extract_greenland_zones.py -- MAINTAINER STEP: reduce the three gridded
temperature products to the annual Greenland-zone and global series that
build_t_gis.py needs, and write them to data/sources/greenland_zones/.

PORT of the measurement half of the research `python/build_t_gis.py`
(SLR-RFF-BRICK @ d59ce12): region_paths, berkeley_land_lookup, cell_weights,
the three loaders, annual_weighted, global_annual, rebase and build() are
copied verbatim. What changed is I/O only:
  * paths come from --raw-dir (the gridded products, the GTN-G shapefile)
    instead of the research tree, and each file is checked against
    tools/targets/inputs.sha256 before it is read (a mismatch aborts);
  * instead of the driver and the amplification tables, it writes the series
    those are computed from, one file per product, at full float64
    round-trip precision (pandas' default float output), so build_t_gis.py
    recomputes every downstream number from exactly the same doubles;
  * the number of grid cells each zone touches (n_cells, a column of the
    research gis_driver_constants.csv) is carried along with the series.

Rows: one per (zone, year) the research build() produced, zone in
south / all / central / north / global. A year missing from a series (fewer
than 12 months, or no weighted cell) is absent, exactly as in the research
pandas Series, so rolling means downstream see the same positions.

The standard rebuild does not run this: it reads the shipped extraction.
Gate: tools/targets/check_targets_identity.py --extraction compares a fresh
extraction with the shipped files byte for byte.

Usage: python extract_greenland_zones.py --raw-dir RAW --out-dir DIR
"""
import argparse
import os
import sys

import netCDF4
import numpy as np
import pandas as pd
import shapefile as pyshp
from matplotlib.path import Path as MplPath

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import target_inputs as TI  # noqa: E402

# ---- inputs (under --raw-dir) ------------------------------------------------
HAD_REL = "HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc"
BERK_REL = "Land_and_Ocean_LatLong1.nc"
GIST_REL = "gistemp1200_GHCNv4_ERSSTv5.nc"
SHP_REL = "GlacReg_2023/GTN-G_202307_o1regions.shp"
HAD_NC = BERK_NC = GIST_NC = SHP = None          # bound in main()

# ---- zones (verbatim from research build_t_gis.py) ---------------------------
GTNG_REGION = "05"                  # "Greenland Periphery" -- the region polygon
ZONES = {
    "south": (59.0, 70.0),          # HEADLINE: the ablation-dominated south
    "all": (59.0, 84.0),            # pre-registered sensitivity arm
}
DIAG_ZONES = {"central": (70.0, 77.0), "north": (77.0, 84.0)}
SUBGRID = 5                         # SUBGRID^2 point-in-polygon samples per cell
LAND_FRAC_FLOOR = 0.0               # land fraction used as a WEIGHT, not a threshold
MIN_MONTHS = 12                     # calendar-year completeness
BASE = (1850, 1900)                 # driver anomaly baseline (multi-year baseline rule)
BASE_COMMON = (1880, 1900)          # baseline all three products can supply (GISTEMP@1880)
PRODUCT_BASELINE = {                # first year each product can be rebased from
    "HadCRUT5": BASE, "BerkeleyEarth": BASE, "GISTEMP": BASE_COMMON,
}
GLOBAL_LABEL = "global"
OUT_NAME = "greenland_zone_series_{product}.csv"


# =============================================================================
# 1. the Greenland-land mask, shared by every product (verbatim)
# =============================================================================
def region_paths(region):
    """Closed rings of one GTN-G first-order region as matplotlib Paths."""
    sf = pyshp.Reader(SHP)
    fields = [f[0] for f in sf.fields[1:]]
    key = fields.index("o1region")
    out = []
    for sr in sf.shapeRecords():
        if str(sr.record[key]).zfill(2) != region:
            continue
        parts = list(sr.shape.parts) + [len(sr.shape.points)]
        for i in range(len(sr.shape.parts)):
            ring = np.asarray(sr.shape.points[parts[i]:parts[i + 1]])
            if len(ring) >= 3:
                out.append(MplPath(ring))
    if not out:
        raise ValueError(f"GTN-G region {region} not found in {SHP}")
    return out


def berkeley_land_lookup():
    """Nearest-neighbour land-fraction lookup on Berkeley's 1-deg mask."""
    d = netCDF4.Dataset(BERK_NC)
    blat = np.asarray(d.variables["latitude"][:], dtype=float)
    blon = np.asarray(d.variables["longitude"][:], dtype=float)
    lm = np.asarray(d.variables["land_mask"][:], dtype=float)

    def land_at(lon, lat):
        i = np.clip(np.searchsorted(blat, lat) - 0, 0, len(blat) - 1)
        i = np.clip(np.rint(lat - blat[0]).astype(int), 0, len(blat) - 1)
        j = np.clip(np.rint(lon - blon[0]).astype(int), 0, len(blon) - 1)
        return lm[i, j]

    return land_at


def cell_weights(lat, lon, rings, land_at, lat_lo, lat_hi):
    """
    Per-cell weight = cos(lat) * (fraction of the cell that is Greenland land
    inside [lat_lo, lat_hi]), by SUBGRID^2 point sampling.
    """
    dlat = float(np.abs(np.diff(lat)).mean())
    dlon = float(np.abs(np.diff(lon)).mean())
    off = (np.arange(SUBGRID) + 0.5) / SUBGRID - 0.5
    w = np.zeros((len(lat), len(lon)))
    for i, la in enumerate(lat):
        sub_la = la + off * dlat
        if sub_la.max() < lat_lo - dlat or sub_la.min() > lat_hi + dlat:
            continue
        for j, lo in enumerate(lon):
            lo180 = ((lo + 180.0) % 360.0) - 180.0
            sub_lo = lo180 + off * dlon
            pts = np.array([(x, y) for y in sub_la for x in sub_lo])
            keep = (pts[:, 1] >= lat_lo) & (pts[:, 1] <= lat_hi)
            if not keep.any():
                continue
            inside = np.zeros(len(pts), dtype=bool)
            for ring in rings:
                inside |= ring.contains_points(pts)
            frac_land = np.array([land_at(x, y) for x, y in pts])
            hit = keep & inside & (frac_land > LAND_FRAC_FLOOR)
            if hit.any():
                w[i, j] = float((frac_land[hit]).sum() / len(pts)) * \
                    float(np.cos(np.deg2rad(la)))
    return w


# =============================================================================
# 2. product loaders -> (values[time, lat, lon], year, month, lat, lon) (verbatim)
# =============================================================================
def load_hadcrut5():
    d = netCDF4.Dataset(HAD_NC)
    lat = np.asarray(d.variables["latitude"][:], dtype=float)
    lon = np.asarray(d.variables["longitude"][:], dtype=float)
    tv = d.variables["time"]
    dates = netCDF4.num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"))
    yr = np.array([x.year for x in dates])
    mo = np.array([x.month for x in dates])
    v = np.ma.filled(d.variables["tas_mean"][:].astype(float), np.nan)
    return v, yr, mo, lat, lon


def load_berkeley():
    d = netCDF4.Dataset(BERK_NC)
    lat = np.asarray(d.variables["latitude"][:], dtype=float)
    lon = np.asarray(d.variables["longitude"][:], dtype=float)
    t = np.asarray(d.variables["time"][:], dtype=float)
    yr = np.floor(t + 1e-6).astype(int)
    mo = np.rint((t - yr) * 12 - 0.5).astype(int) + 1
    v = np.ma.filled(d.variables["temperature"][:].astype(float), np.nan)
    return v, yr, mo, lat, lon


def load_gistemp():
    d = netCDF4.Dataset(GIST_NC)
    lat = np.asarray(d.variables["lat"][:], dtype=float)
    lon = np.asarray(d.variables["lon"][:], dtype=float)
    tv = d.variables["time"]
    dates = netCDF4.num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"))
    yr = np.array([x.year for x in dates])
    mo = np.array([x.month for x in dates])
    v = np.ma.filled(d.variables["tempanomaly"][:].astype(float), np.nan)
    return v, yr, mo, lat, lon


PRODUCTS = {"HadCRUT5": load_hadcrut5, "BerkeleyEarth": load_berkeley,
            "GISTEMP": load_gistemp}


# =============================================================================
# 3. weighted annual series (verbatim)
# =============================================================================
def annual_weighted(v, yr, mo, w):
    """Calendar-year mean of the w-weighted spatial mean; NaN cells drop out of
    both numerator and denominator, so partial coverage is handled correctly."""
    wf = w[None, :, :]
    ok = np.isfinite(v)
    num = np.nansum(np.where(ok, v, 0.0) * wf, axis=(1, 2))
    den = np.sum(np.where(ok, 1.0, 0.0) * wf, axis=(1, 2))
    m = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    df = pd.DataFrame({"year": yr, "month": mo, "v": m})
    g = df.dropna().groupby("year")["v"]
    s = g.mean()
    return s[g.count() >= MIN_MONTHS]


def global_annual(v, yr, mo, lat, lon):
    w = np.repeat(np.cos(np.deg2rad(lat))[:, None], len(lon), axis=1)
    return annual_weighted(v, yr, mo, w)


def rebase(s, base):
    seg = s.loc[base[0]:base[1]]
    if len(seg) == 0:
        raise ValueError("empty baseline segment")
    return s - seg.mean()


# =============================================================================
# 4. build (verbatim)
# =============================================================================
def build():
    rings = region_paths(GTNG_REGION)
    for name, (x, y), want in [("Iceland", (-19.0, 65.0), False),
                               ("Baffin", (-70.0, 68.0), False),
                               ("Ellesmere", (-75.0, 79.0), False),
                               ("S-Greenland", (-45.0, 62.0), True),
                               ("C-Greenland", (-40.0, 72.0), True)]:
        got = any(r.contains_point((x, y)) for r in rings)
        assert got == want, f"GTN-G region {GTNG_REGION} mask: {name} inside={got}"
    land_at = berkeley_land_lookup()

    series, globals_, weights = {}, {}, {}
    all_zones = {**ZONES, **DIAG_ZONES}
    for pname, loader in PRODUCTS.items():
        v, yr, mo, lat, lon = loader()
        base = PRODUCT_BASELINE[pname]
        globals_[pname] = rebase(global_annual(v, yr, mo, lat, lon), base)
        for zname, (lo, hi) in all_zones.items():
            w = cell_weights(lat, lon, rings, land_at, lo, hi)
            weights[(pname, zname)] = w
            series[(pname, zname)] = rebase(annual_weighted(v, yr, mo, w), base)
        print(f"  {pname:14s} grid {len(lat)}x{len(lon)}  "
              + "  ".join(f"{z}: {int((weights[(pname, z)] > 0).sum())} cells"
                          for z in ZONES), flush=True)
    return series, globals_, weights


# =============================================================================
# 5. write (new)
# =============================================================================
def to_long(pname, series, globals_, weights):
    parts = []
    for zname in list(ZONES) + list(DIAG_ZONES):
        s = series[(pname, zname)]
        parts.append(pd.DataFrame({"zone": zname, "n_cells": int((weights[(pname, zname)] > 0).sum()),
                                   "year": s.index.astype(int), "value": s.to_numpy(dtype=float)}))
    g = globals_[pname]
    parts.append(pd.DataFrame({"zone": GLOBAL_LABEL, "n_cells": 0,
                               "year": g.index.astype(int), "value": g.to_numpy(dtype=float)}))
    return pd.concat(parts, ignore_index=True)


def main():
    global HAD_NC, BERK_NC, GIST_NC, SHP
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raw-dir", required=True, help="directory holding the raw products (see inputs.sha256)")
    ap.add_argument("--out-dir", required=True, help="where to write greenland_zone_series_<product>.csv")
    ap.add_argument("--products", nargs="+", choices=list(PRODUCTS), default=list(PRODUCTS),
                    help="products to extract (default all three); the Berkeley land fraction is always read")
    ap.add_argument("--no-verify", action="store_true",
                    help="skip the sha256 check of the raw inputs (the output then cannot be compared)")
    a = ap.parse_args()
    HAD_NC, BERK_NC, GIST_NC, SHP = (os.path.join(a.raw_dir, r) for r in (HAD_REL, BERK_REL, GIST_REL, SHP_REL))
    for p in [q for q in PRODUCTS if q not in a.products]:
        del PRODUCTS[p]
    if not a.no_verify:
        pins = TI.read_pins()
        need = {"HadCRUT5": HAD_REL, "BerkeleyEarth": BERK_REL, "GISTEMP": GIST_REL}
        for rel in sorted({need[p] for p in PRODUCTS} | {BERK_REL}) + [SHP_REL, SHP_REL[:-4] + ".shx",
                                                                       SHP_REL[:-4] + ".dbf"]:
            got = TI.sha256(os.path.join(a.raw_dir, rel))
            if got != pins[f"raw/{rel}"]:
                sys.exit(f"ABORT: {rel} sha256 {got} != pinned {pins['raw/' + rel]} (inputs.sha256)")
        print("raw inputs match inputs.sha256")
    series, globals_, weights = build()
    os.makedirs(a.out_dir, exist_ok=True)
    for pname in PRODUCTS:
        out = os.path.join(a.out_dir, OUT_NAME.format(product=pname))
        to_long(pname, series, globals_, weights).to_csv(out, index=False)   # default float output: round-trip
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
