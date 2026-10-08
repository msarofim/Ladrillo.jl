"""Identity gate 5 (Table 5): python/ic_ladrillo_vs_brick20.py against the research outputs, for each rho bound.

    python test/identity/compare_ic_table.py <research outputs dir> [rho ...]     (default rho: 0.99 0.95 0.9)

Compares ic_ladrillo_vs_brick20_L27{sfx}.csv with outputs/tables/ic_ladrillo_vs_brick20{sfx}.csv on every
(arm, model, series, stat) value, and the _perdraw files on every column, with exact equality on the parsed floats
(NaN == NaN), parsed with float_precision="round_trip" (pandas' default parser is not exact and can merge values
1 ulp apart). `provenance` is excluded: it names the producing script and its inputs, which differ by construction.
"""
import os
import sys

import numpy as np
import pandas as pd

RES = os.path.abspath(sys.argv[1])
RHOS = [float(r) for r in sys.argv[2:]] or [0.99, 0.95, 0.9]
PKG = os.environ.get("LADRILLO_TABLES") or os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "outputs", "tables")
KEY = ["arm", "model", "series", "stat"]

bad = 0
for rho in RHOS:
    sfx = "" if rho == 0.99 else f"_rho{rho:g}"
    for kind in ("", "_perdraw"):
        a = pd.read_csv(os.path.join(RES, f"ic_ladrillo_vs_brick20_L27{sfx}{kind}.csv"), float_precision="round_trip")
        b = pd.read_csv(os.path.join(PKG, f"ic_ladrillo_vs_brick20{sfx}{kind}.csv"), float_precision="round_trip")
        a, b = (d.drop(columns=["provenance"], errors="ignore") for d in (a, b))
        name = f"rho {rho:<5g} {'perdraw' if kind else 'table'}"
        if list(a.columns) != list(b.columns) or len(a) != len(b):
            print(f"  {name:22s} SHAPE: {a.shape} vs {b.shape}, cols equal {list(a.columns) == list(b.columns)}")
            bad += 1
            continue
        if not kind and not (a[KEY].astype(str).values == b[KEY].astype(str).values).all():
            print(f"  {name:22s} KEY ORDER differs")
            bad += 1
            continue
        num = [c for c in a.columns if c not in KEY and c != "model"]
        x, y = a[num].to_numpy(float), b[num].to_numpy(float)
        diff = ~((x == y) | (np.isnan(x) & np.isnan(y)))
        nd = int(diff.sum())
        w = float(np.nanmax(np.abs(x - y)[diff])) if nd else 0.0
        if kind and not (a["model"].values == b["model"].values).all():
            nd += 1
        bad += nd
        print(f"  {name:22s} {x.size:9d} values  {nd} differ  max|d| {w:.3e}  {'IDENTICAL' if nd == 0 else 'FAIL'}")
print(f"[IC-TABLE-IDENTITY] {'PASS' if bad == 0 else 'FAIL'}")
sys.exit(0 if bad == 0 else 1)
