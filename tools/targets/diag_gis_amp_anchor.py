#!/usr/bin/env python3
"""
diag_gis_amp_anchor.py -- at what warming level does the observed Greenland
amplification sit? (the anchor of the amplification shape S(dT))

PORT of the research `python/diag_gis_amp_anchor.py` (SLR-RFF-BRICK @ 8b68d19).
amp_and_anchor() and the table are verbatim. What changed (I/O only): the
series come from build_t_gis.load() (the shipped extraction) instead of
re-running the research build_t_gis.build() on the gridded products, paths
come from --sources / --work, and the headline printout is dropped.

The observed amplification is a through-origin fit over 1901-2024,
amp = sum(x*y)/sum(x^2), x the product's global anomaly and y the south-
Greenland anomaly, both rel 1850-1900. Written as a weighted mean of y/x the
weights are x^2, so the warming level it represents is
dT_eff = sum(x^3)/sum(x^2). That is where S must equal 1.

Gate (kept from the research script): the recomputed amp must reproduce
outputs/gis_driver_constants.csv (written by build_t_gis.py in the same
build) to 1e-9.

Reads:  --sources greenland_zones/ ; --work outputs/gis_driver_constants.csv
Writes: --work outputs/gis_amp_anchor.csv   per product x zone x window: amp, dT_eff, n

Usage: python diag_gis_amp_anchor.py --sources DIR --work DIR
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_t_gis as BTG  # noqa: E402

OUT_REL = "outputs/gis_amp_anchor.csv"
AMP_TOL = 1e-9                             # gate: recomputed amp vs committed


def amp_and_anchor(zone, glob, win):
    """(amp, dT_eff, n) for one product x zone x window."""
    x = glob.loc[win[0]:win[1]]
    y = zone.reindex(x.index)
    ok = x.notna() & y.notna()
    if ok.sum() < 20:
        return np.nan, np.nan, int(ok.sum())
    xv, yv = x[ok].to_numpy(float), y[ok].to_numpy(float)
    sxx = (xv ** 2).sum()
    return float((xv * yv).sum() / sxx), float((xv ** 3).sum() / sxx), int(ok.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sources", required=True)
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    OUT_CSV = os.path.join(a.work, OUT_REL)
    CONSTANTS_CSV = os.path.join(a.work, BTG.OUT_CONSTANTS_REL)

    series, globals_, _ = BTG.load(a.sources)
    committed = pd.read_csv(CONSTANTS_CSV)

    rows = []
    for pname in BTG.PRODUCTS:
        for zname in list(BTG.ZONES) + list(BTG.DIAG_ZONES):
            for wname, win in BTG.AMP_WINDOWS.items():
                a_, dte, n = amp_and_anchor(series[(pname, zname)], globals_[pname], win)
                x = globals_[pname].loc[win[0]:win[1]]
                y = series[(pname, zname)].reindex(x.index)
                ok = x.notna() & y.notna()
                rows.append(dict(product=pname, zone=zname, window=wname,
                                 year0=win[0], year1=win[1], n_years=n,
                                 amp=a_, dt_eff=dte,
                                 dt_mean=float(x[ok].mean()),
                                 dt_last=float(x[ok].iloc[-1])))
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)

    bad = []
    for _, r in df.iterrows():
        c = committed[(committed["product"] == r["product"])
                      & (committed.zone == r.zone)]
        if c.empty:
            continue
        want = float(c[f"amp_{r.window}"].iloc[0])
        if not (np.isnan(want) and np.isnan(r.amp)) and abs(want - r.amp) > AMP_TOL:
            bad.append(f"{r['product']}/{r.zone}/{r.window}: {r.amp:.12f} vs {want:.12f}")
    if bad:
        sys.exit("GATE FAIL -- recomputed amp does not match "
                 f"{os.path.basename(CONSTANTS_CSV)}:\n  " + "\n  ".join(bad))
    h = df[(df.zone == BTG.HEADLINE_ZONE) & (df.window == BTG.AMP_WINDOW_HEADLINE)]
    print(f"diag_gis_amp_anchor | amp reproduces gis_driver_constants.csv to {AMP_TOL:g} "
          f"({len(df)} cells) | anchor dT_eff {BTG.HEADLINE_ZONE}/{BTG.AMP_WINDOW_HEADLINE} "
          f"= {h.dt_eff.mean():.4f} K | wrote {OUT_REL}")


if __name__ == "__main__":
    main()
