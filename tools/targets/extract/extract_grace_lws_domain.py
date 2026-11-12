#!/usr/bin/env python3
"""
extract_grace_lws_domain.py -- MAINTAINER STEP: the monthly GRACE/GRACE-FO mass
anomaly over global land minus both ice sheets, i.e.
data/sources/grace/grace_lws_domain_monthly.csv.

PORT of the gridded half of the research `python/build_lws_grace_extension.py`
(SLR-RFF-BRICK @ d2e23d3): the grid areas, region_mask, mass_gt, the ice-sheet
masks and the five gates ([AREA], [ANT-SWEEP], [GRL]/[ANT], [DOMAIN], [SCOPE])
are verbatim. What changed (I/O only):
  * the mascon file is named, not found by glob. The research script took the
    NEWEST file matching GRCTellus.JPL.*.GLO.RL06.3M.MSCNv04CRI.nc, so its
    result changed whenever a newer month was downloaded next to it. The file
    here is the one inputs.sha256 pins (2002-04 to 2026-06);
  * paths from --raw-dir (mascon, GTN-G shapefile) and --sources (the two JPL
    ice-sheet series the masks are gated against), all checked against
    inputs.sha256;
  * it writes the monthly domain mass at full round-trip precision and stops.
    Deseasonalising, annual means, the glacier subtraction and the splice are
    in tools/targets/build_lws_grace_extension.py.

Domain: land cells (the mascon land_mask) minus GTN-G region 05 (Greenland ice
sheet + periphery) minus Antarctica (GTN-G region 20 plus every land cell south
of 60 S, because the region-20 polygon misses some cells at 0.5 degree).
Units: Gt (1 cm water over 1 km2 = 1e-5 Gt).

Usage: python extract_grace_lws_domain.py --raw-dir RAW --sources data/sources --out-dir DIR
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import shapefile
import xarray as xr
from matplotlib.path import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import target_inputs as TI  # noqa: E402

MASCON_REL = "GRCTellus.JPL.200204_202606.GLO.RL06.3M.MSCNv04CRI.nc"
GTNG_REL = "GlacReg_2023/GTN-G_202307_o1regions"
OUT_NAME = "grace_lws_domain_monthly.csv"

# ---- physical / grid constants (verbatim) --------------------------------------
R_EARTH_KM = 6371.0072
CM_KM2_TO_GT = 1e-5        # 1 cm of water over 1 km2 = 1e4 m3 = 1e7 kg = 1e-5 Gt.
EARTH_AREA_KM2 = 5.101e8   # for the [AREA] gate
# ---- masking / scope (verbatim) ------------------------------------------------
GIS_REGION = "05"          # GTN-G Greenland; the polygon covers the ice sheet + periphery
AIS_REGION = "20"          # GTN-G Antarctic Mainland
AIS_SWEEP_LAT = -60.0      # the region-20 polygon MISSES land cells; see [ANT-SWEEP]
# ---- gate thresholds (verbatim) -----------------------------------------------
ICE_GATE_TOL = 0.02        # |ratio-1| against JPL's own published ice-sheet series
AREA_GATE_TOL = 0.001


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--sources", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    pins = TI.read_pins()
    for rel in (MASCON_REL, GTNG_REL + ".shp", GTNG_REL + ".shx", GTNG_REL + ".dbf"):
        got = TI.sha256(os.path.join(a.raw_dir, rel))
        if got != pins["raw/" + rel]:
            sys.exit(f"ABORT: {rel} sha256 {got} != pinned {pins['raw/' + rel]} (inputs.sha256)")
    for rel in ("grace/grace_greenland_mass.txt", "grace/grace_antarctica_mass.txt"):
        if TI.sha256(os.path.join(a.sources, rel)) != pins["sources/" + rel]:
            sys.exit(f"ABORT: {rel} does not match inputs.sha256")
    GTNG = os.path.join(a.raw_dir, GTNG_REL)

    # ---- verbatim from here (build_lws_grace_extension.py:105-191) ----
    ds = xr.open_dataset(os.path.join(a.raw_dir, MASCON_REL))
    print(f"  file     {MASCON_REL} ({ds.sizes['time']} months)")
    lat, lon = ds.lat.values, ds.lon.values
    lon180 = np.where(lon > 180, lon - 360, lon)
    LON, LAT = np.meshgrid(lon180, lat)
    pts = np.column_stack([LON.ravel(), LAT.ravel()])
    land = ds.land_mask.values == 1
    lwe = ds.lwe_thickness.values
    lb = ds.lat_bounds.values
    dlon_deg = float(abs(np.diff(np.sort(lon))[0]))
    area = (R_EARTH_KM ** 2 * np.deg2rad(dlon_deg)
            * (np.sin(np.deg2rad(lb[:, 1])) - np.sin(np.deg2rad(lb[:, 0]))))[:, None] \
        * np.ones((1, len(lon)))
    t = pd.to_datetime(ds.time.values)
    tdec = t.year + (t.dayofyear - 1) / 365.25

    def region_mask(code):
        """Union of all polygon parts of one GTN-G o1 region, on the mascon grid."""
        sf = shapefile.Reader(GTNG)
        m = np.zeros(len(pts), bool)
        for s, r in zip(sf.shapes(), sf.records()):
            if str(r[0]) != code:
                continue
            parts = list(s.parts) + [len(s.points)]
            for a_, b_ in zip(parts[:-1], parts[1:]):
                if b_ - a_ >= 3:
                    m |= Path(np.asarray(s.points[a_:b_])).contains_points(pts)
        return m.reshape(LAT.shape)

    def mass_gt(mask):
        """Mass anomaly time series over a mask, Gt."""
        return np.nansum(lwe * (area * mask)[None, :, :], axis=(1, 2)) * CM_KM2_TO_GT

    ok_area = abs(area.sum() / EARTH_AREA_KM2 - 1) < AREA_GATE_TOL
    print(f"[AREA]  grid sums to {area.sum():.4e} km2 vs {EARTH_AREA_KM2:.4e} -> "
          f"{'PASS' if ok_area else 'FAIL'}")
    assert ok_area

    m_gis = region_mask(GIS_REGION) & land
    m_ais_poly = region_mask(AIS_REGION) & land
    m_south = land & (lat[:, None] < AIS_SWEEP_LAT)
    m_ais = m_ais_poly | m_south

    missed = m_south & ~m_ais_poly
    A = np.vstack([np.ones(len(tdec)), tdec - tdec[0]]).T
    tr_missed = np.linalg.lstsq(A, mass_gt(missed), rcond=None)[0][1]
    print(f"[ANT-SWEEP] the region-{AIS_REGION} polygon misses {missed.sum()} land cells worth "
          f"{tr_missed:+.2f} Gt/yr -> swept in via lat < {AIS_SWEEP_LAT}")

    for code, name, ref, mask in (("GRL", "Greenland", "grace/grace_greenland_mass.txt", m_gis),
                                  ("ANT", "Antarctica", "grace/grace_antarctica_mass.txt", m_ais)):
        r = np.loadtxt(os.path.join(a.sources, ref), comments="HDR")
        idx = np.array([np.argmin(abs(r[:, 0] - x)) for x in tdec])
        sel = np.array([abs(r[i, 0] - x) < 0.02 for i, x in zip(idx, tdec)])
        v = mass_gt(mask)
        ours = v[sel] - v[sel][0]
        theirs = r[idx[sel], 1] - r[idx[sel], 1][0]
        Ag = np.vstack([np.ones(sel.sum()), tdec[sel] - tdec[sel][0]]).T
        bo = np.linalg.lstsq(Ag, ours, rcond=None)[0][1]
        bt = np.linalg.lstsq(Ag, theirs, rcond=None)[0][1]
        good = abs(bo / bt - 1) < ICE_GATE_TOL
        print(f"[{code}]   {name:10s} our trend {bo:8.2f} vs JPL {bt:8.2f} Gt/yr | ratio "
              f"{bo / bt:.4f} | rms {np.sqrt(((ours - theirs) ** 2).mean()):5.1f} Gt -> "
              f"{'PASS' if good else 'FAIL'}")
        assert good, f"{name} mask does not reproduce the published series"

    ice = m_gis | m_ais
    dom = land & ~ice
    print(f"[DOMAIN] LWS domain {dom.sum()} cells, {area[dom].sum() / 1e6:.2f} Mkm2 "
          f"(global land minus both ice sheets)")

    n19 = (region_mask("19") & dom).sum()
    print(f"[SCOPE] region-19 land cells inside the domain: {n19} -> "
          f"{'PASS (nothing to subtract)' if n19 == 0 else 'FAIL'}")
    assert n19 == 0

    g = mass_gt(dom)
    # ---- end verbatim ----

    stamps = [x.isoformat() for x in t]
    assert pd.DatetimeIndex(pd.to_datetime(stamps)).equals(t), "time stamps do not round-trip"
    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, OUT_NAME)
    pd.DataFrame({"time": stamps, "domain_mass_gt": g}).to_csv(out, index=False)   # round-trip floats
    print(f"extract_grace_lws_domain | {len(t)} months | wrote {out}")


if __name__ == "__main__":
    main()
