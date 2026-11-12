#!/usr/bin/env python3
"""
build_t_gis.py -- the Greenland ice-sheet temperature driver for Ladrillo and
its amplification prior.

PORT of the derivation half of the research `python/build_t_gis.py`
(SLR-RFF-BRICK @ d59ce12). The research script measured the zone and global
series on the gridded products and then derived everything below from them.
Here the measurement half lives in extract/extract_greenland_zones.py (a
maintainer step that needs the gridded products), and this script starts
from its shipped output, data/sources/greenland_zones/. Everything from the
series on is copied verbatim: amp_through_origin, melt_rate_target,
corr_with_meltrate, amp_prior and main()'s tables.

What changed (I/O only):
  * paths from --sources / --work;
  * the series come from the extraction (read at round-trip precision) instead
    of the gridded products;
  * the DRIVER GUARD is gone. The research script read its own previous output
    and refused to overwrite a changed driver; that is a gate reading the path
    its script writes. The byte-identity gate (check_targets_identity.py)
    compares the rebuilt driver with the shipped one instead;
  * t_gis_zones_allproducts.csv, the provenance markdown, the console tables
    and the figure are dropped: nothing downstream reads them.

Construction (unchanged; see the research docstring for the reasoning):
  zone south = 59-70 N (headline), all = 59-84 N, central = 70-77 N,
  north = 77-84 N; mask = GTN-G 2023 region 05 x Berkeley Earth 1-degree land
  fraction; annual = calendar-year mean of 12 months; anomalies rel 1850-1900
  (GISTEMP 1880-1900); amplification = through-origin fit of the zone on the
  same product's global mean, windows full 1901-2024, early 1901-1960, modern
  1961-2024; prior = mean and population sd over the three products.

Reads (under --sources): greenland_zones/greenland_zone_series_<product>.csv
       (under --work):   outputs/recalib_targets_ext.csv (the GIS target, for the
                         melt-rate correlation column only)
Writes (under --work):
  data/observations/t_gis_zones.csv   the driver, HadCRUT5, one column per zone (%.12f)
  outputs/gis_driver_constants.csv    amplification + melt-rate correlation per product x zone
  outputs/gis_amp_prior.csv           cross-product amplification prior per zone x window

Usage: python build_t_gis.py --sources DIR --work DIR
"""
import argparse
import os

import numpy as np
import pandas as pd

# ---- zones -----------------------------------------------------------------
ZONES = {
    "south": (59.0, 70.0),          # HEADLINE: the ablation-dominated south
    "all": (59.0, 84.0),            # pre-registered sensitivity arm
}
DIAG_ZONES = {"central": (70.0, 77.0), "north": (77.0, 84.0)}
DRIVER_ZONES = {**ZONES, **DIAG_ZONES}      # columns written to t_gis_zones.csv
HEADLINE_ZONE = "south"

# ---- construction constants ------------------------------------------------
BASE = (1850, 1900)                 # driver anomaly baseline (multi-year baseline rule)
BASE_COMMON = (1880, 1900)          # baseline all three products can supply (GISTEMP@1880)
AMP_WINDOWS = {"full": (1901, 2024), "early": (1901, 1960), "modern": (1961, 2024)}
AMP_WINDOW_HEADLINE = "full"        # matches ladrillo_data.AMP_FIT_WIN
CORR_WIN = (1900, 2018)             # the Frederikse GIS target window
SMOOTH = 11                         # melt-rate smoothing, years

PRODUCT_BASELINE = {                # first year each product can be rebased from
    "HadCRUT5": BASE, "BerkeleyEarth": BASE, "GISTEMP": BASE_COMMON,
}
PRODUCTS = list(PRODUCT_BASELINE)   # the research PRODUCTS dict order
GLOBAL_LABEL = "global"
SERIES_REL = "greenland_zones/greenland_zone_series_{product}.csv"
TARGETS_REL = "outputs/recalib_targets_ext.csv"
OUT_DRIVER_REL = "data/observations/t_gis_zones.csv"
OUT_CONSTANTS_REL = "outputs/gis_driver_constants.csv"
OUT_PRIORS_REL = "outputs/gis_amp_prior.csv"
DRIVER_PRODUCT = "HadCRUT5"
DRIVER_FLOAT_FORMAT = "%.12f"

TARGETS_CSV = None                  # bound in main() / load()


# =============================================================================
# the series, from the shipped extraction
# =============================================================================
def load(sources):
    """(series, globals_, n_cells) exactly as the research build() returned them:
    series[(product, zone)] and globals_[product] are pandas Series indexed by year."""
    series, globals_, n_cells = {}, {}, {}
    for pname in PRODUCTS:
        df = pd.read_csv(os.path.join(sources, SERIES_REL.format(product=pname)),
                         float_precision="round_trip")
        for zname, g in df.groupby("zone", sort=False):
            s = pd.Series(g["value"].to_numpy(dtype=float), index=pd.Index(g["year"].to_numpy(), name="year"),
                          name="v")
            if zname == GLOBAL_LABEL:
                globals_[pname] = s
            else:
                series[(pname, zname)] = s
                n_cells[(pname, zname)] = int(g["n_cells"].iloc[0])
    missing = [(p, z) for p in PRODUCTS for z in DRIVER_ZONES if (p, z) not in series]
    assert not missing and set(globals_) == set(PRODUCTS), f"extraction incomplete: {missing}"
    return series, globals_, n_cells


# =============================================================================
# 3. (verbatim) amplification and relevance
# =============================================================================
def amp_through_origin(zone, glob, win):
    x = glob.loc[win[0]:win[1]]
    y = zone.reindex(x.index)
    ok = x.notna() & y.notna()
    if ok.sum() < 20:
        return np.nan
    return float((x[ok] * y[ok]).sum() / (x[ok] ** 2).sum())


def melt_rate_target():
    tgt = pd.read_csv(TARGETS_CSV).set_index("year")
    return tgt["gis"].dropna().diff().rolling(SMOOTH, center=True).mean()


def corr_with_meltrate(s, melt_rate):
    win = np.arange(CORR_WIN[0], CORR_WIN[1] + 1)
    sm = s.rolling(SMOOTH, center=True).mean().reindex(win)
    both = pd.concat([sm, melt_rate.reindex(win)], axis=1).dropna()
    if len(both) <= 30:
        return np.nan
    return float(np.corrcoef(both.iloc[:, 0], both.iloc[:, 1])[0, 1])


def amp_prior(con, zone, window):
    """Cross-product amplification prior for one zone: mean, sd, and the
    product range it spans. sd is the population sd over the three products --
    the products are the uncertainty here, and there are only three of them."""
    v = con[con.zone == zone][f"amp_{window}"].to_numpy(dtype=float)
    return dict(zone=zone, window=window, mean=float(v.mean()),
                sd=float(v.std(ddof=0)), lo=float(v.min()), hi=float(v.max()),
                spread_ratio=float(v.max() / v.min()))


def main():
    global TARGETS_CSV
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sources", required=True)
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    TARGETS_CSV = os.path.join(a.work, TARGETS_REL)
    OUT_DRIVER = os.path.join(a.work, OUT_DRIVER_REL)
    OUT_CONSTANTS = os.path.join(a.work, OUT_CONSTANTS_REL)
    OUT_PRIORS = os.path.join(a.work, OUT_PRIORS_REL)
    for p in (OUT_DRIVER, OUT_CONSTANTS, OUT_PRIORS):
        os.makedirs(os.path.dirname(p), exist_ok=True)

    series, globals_, n_cells = load(a.sources)
    melt_rate = melt_rate_target()

    rows = []
    for (pname, zname), s in series.items():
        r = dict(product=pname, zone=zname,
                 baseline=f"{PRODUCT_BASELINE[pname][0]}-{PRODUCT_BASELINE[pname][1]}",
                 year_first=int(s.index.min()), year_last=int(s.index.max()),
                 n_cells=n_cells[(pname, zname)],
                 corr_meltrate=corr_with_meltrate(s, melt_rate))
        for wname, win in AMP_WINDOWS.items():
            r[f"amp_{wname}"] = amp_through_origin(s, globals_[pname], win)
        rows.append(r)
    con = pd.DataFrame(rows).sort_values(["zone", "product"]).reset_index(drop=True)
    con.to_csv(OUT_CONSTANTS, index=False)

    priors = pd.DataFrame([amp_prior(con, z, w)
                           for z in list(ZONES) + list(DIAG_ZONES)
                           for w in AMP_WINDOWS])
    priors.to_csv(OUT_PRIORS, index=False)

    # headline driver: every zone, one file, full precision
    drv = pd.concat({z: series[(DRIVER_PRODUCT, z)] for z in DRIVER_ZONES}, axis=1)
    drv.index.name = "year"
    drv.to_csv(OUT_DRIVER, float_format=DRIVER_FLOAT_FORMAT)

    pr = amp_prior(con, HEADLINE_ZONE, AMP_WINDOW_HEADLINE)
    print(f"build_t_gis | driver {DRIVER_PRODUCT} {int(drv.index.min())}-{int(drv.index.max())} | "
          f"amp prior {HEADLINE_ZONE}/{AMP_WINDOW_HEADLINE} N({pr['mean']:.3f}, {pr['sd']:.3f})")
    for p in (OUT_DRIVER_REL, OUT_CONSTANTS_REL, OUT_PRIORS_REL):
        print(f"  wrote {p}")


if __name__ == "__main__":
    main()
