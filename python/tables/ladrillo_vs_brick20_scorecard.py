#!/usr/bin/env python3
"""
ladrillo_vs_brick20_scorecard.py -- Table 4: Ladrillo vs BRICK 2.0 on the hindcast, component by component, in one
set of metrics (mean bias, RMSE, largest absolute error, 90 % band coverage) and in five windows.

PORT of SLR-RFF-BRICK python/scope_ladrillo_vs_brick20_scorecard.py (the research file is
outputs/scope_ladrillo_vs_brick20_scorecard_L27.csv). The windows, components, metrics and the re-reference check are
copied unchanged; what differs is the I/O: one posterior (L27), so there is no --tag, and the inputs are this
repository's own products.

  python python/tables/ladrillo_vs_brick20_scorecard.py [--outdir=outputs/tables]
         [--ladrillo=outputs/hindcast/postpred_components_timeseries.csv]
         [--brick20=outputs/hindcast/postpred_brick20_components_timeseries.csv]
         [--targets=data/observations/calibration_targets.csv]

INPUTS: the two posterior-predictive band files (scripts/posterior_predictive.jl and
scripts/posterior_predictive_brick20.jl) and the calibration targets. The input options exist so that the test suite
can run this script on the research references in test/reference/.

BOTH ARMS ARE ON THE SAME FOOTING (checked, not assumed):
  * the same re-reference window, 1995-2005, and the targets are zeroed on it (the script stops if they are not), so
    the model bands and the observations share one baseline;
  * the same FaIR mean GMST and ocean heat forcing;
  * both bands are PARAMETRIC (no AR(1) noise added), so `coverage90` is comparable to the Ladrillo posterior
    predictive's parametric coverage, not to its predictive coverage.

BRICK 2.0 here is MimiBRICK v2.0.0 on its published joint posterior (Zenodo 20592337), as its authors calibrated it.
It was not recalibrated on these targets, so part of its bias against them is target vintage; that favours Ladrillo,
most of all in the earliest, sparsest window.

Writes <outdir>/ladrillo_vs_brick20_scorecard.csv: one row per (window, component, arm), plus a provenance column.
python/tables/check_tables_identity.py gates it against the research file.
"""
import argparse
import os
import subprocess
from datetime import datetime

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = os.path.relpath(os.path.abspath(__file__), REPO)
TAG = "L27"                                        # the posterior's calibration (every version)
DEFAULT_OUTDIR = os.path.join(REPO, "outputs", "tables")
DEFAULT_LADRILLO = os.path.join(REPO, "outputs", "hindcast", "postpred_components_timeseries.csv")
DEFAULT_BRICK20 = os.path.join(REPO, "outputs", "hindcast", "postpred_brick20_components_timeseries.csv")
DEFAULT_TARGETS = os.path.join(REPO, "data", "observations", "calibration_targets.csv")
OUT_STEM = "ladrillo_vs_brick20_scorecard"

REF = (1995, 2005)                                 # the shared re-reference window
ARMS = ("Ladrillo", "BRICK 2.0")
## component label -> (Ladrillo stem, BRICK 2.0 stem, target column; None = use the Ladrillo file's <stem>_obs)
COMPONENTS = [("AIS", "ais", "ais", "ais"),
              ("glaciers", "glaciers", "gsic", "gsic"),
              ("Greenland", "gis", "gis", "gis"),
              ("thermal exp.", "te", "te", "steric"),
              ("TOTAL", "total", "total", None)]
## Lower band column per arm: the two drivers name their 5th percentile differently.
LOWER = {"Ladrillo": "p05", "BRICK 2.0": "p5"}
## Windows. 1950 and 1993 separate the sparse-observation era, the pre-satellite era and the altimetry era; "full" is
## every year both arms cover.
WINDOWS = [("full", None), ("1900-1919", (1900, 1919)), ("1920-1949", (1920, 1949)),
           ("1950-1992", (1950, 1992)), ("1993-2025", (1993, 2025))]
UNITS = f"cm rel {REF[0]}-{REF[1]}"


def load(lad_path, b20_path, tgt_path):
    a = pd.read_csv(lad_path).set_index("year")
    b = pd.read_csv(b20_path).set_index("year")
    t = pd.read_csv(tgt_path).set_index("year")
    base = t.loc[REF[0]:REF[1], ["ais", "gsic", "gis", "steric"]].mean()
    if base.abs().max() > 1e-3:
        raise SystemExit(f"targets are not zeroed on the shared re-reference window {REF}: {base.to_dict()}; "
                         "the two arms would not be on one baseline")
    return a, b, t


def score(model_p50, lo, hi, obs, window):
    m = obs.notna() & model_p50.notna()
    if window is not None:
        m &= (obs.index >= window[0]) & (obs.index <= window[1])
    if m.sum() == 0:
        return None
    r = (model_p50[m] - obs[m])
    inband = ((obs[m] >= lo[m]) & (obs[m] <= hi[m])).mean()
    return dict(n=int(m.sum()), mean_bias=float(r.mean()),
                rmse=float(np.sqrt((r ** 2).mean())),
                max_abs=float(r.abs().max()), coverage90=float(inband))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    ap.add_argument("--ladrillo", default=DEFAULT_LADRILLO)
    ap.add_argument("--brick20", default=DEFAULT_BRICK20)
    ap.add_argument("--targets", default=DEFAULT_TARGETS)
    args = ap.parse_args()
    outdir = os.path.abspath(args.outdir)

    a, b, t = load(args.ladrillo, args.brick20, args.targets)
    years = a.index.intersection(b.index)
    print("Table 4: Ladrillo vs BRICK 2.0, hindcast scorecard")
    print(f"  overlap {years.min()}-{years.max()}, re-reference {REF[0]}-{REF[1]}, parametric 90% bands both arms\n")

    frames = {"Ladrillo": a, "BRICK 2.0": b}
    rows = []
    for wname, win in WINDOWS:
        print(f"  === window: {wname} ===")
        print(f"  {'component':13s} {'arm':10s} {'n':>4s} {'mean bias':>10s} "
              f"{'RMSE':>8s} {'max|err|':>9s} {'90% cov':>8s}")
        for label, la, ba, tcol in COMPONENTS:
            obs_l = a[f"{la}_obs"] if f"{la}_obs" in a else t[tcol].reindex(a.index)
            obs = obs_l.reindex(years)
            for arm in ARMS:
                df, stem = frames[arm], (la if arm == "Ladrillo" else ba)
                s = score(df[f"{stem}_p50"].reindex(years), df[f"{stem}_{LOWER[arm]}"].reindex(years),
                          df[f"{stem}_p95"].reindex(years), obs, win)
                if s is None:
                    continue
                rows.append(dict(window=wname, component=label, arm=arm, **s))
                print(f"  {label:13s} {arm:10s} {s['n']:4d} {s['mean_bias']:+10.3f} {s['rmse']:8.3f} "
                      f"{s['max_abs']:9.3f} {s['coverage90']:8.1%}")
        print()

    out = pd.DataFrame(rows)
    commit = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    out["provenance"] = (f"Ladrillo {SCRIPT} | posterior {TAG} | seed n/a (no RNG in this script) | "
                         f"ladrillo {os.path.basename(args.ladrillo)} | brick20 {os.path.basename(args.brick20)} | "
                         f"targets {os.path.basename(args.targets)} | commit {commit} | "
                         f"{datetime.now():%Y-%m-%d %H:%M} | {UNITS}; RMSE ratio = Ladrillo/BRICK2.0")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{OUT_STEM}.csv")
    out.to_csv(path, index=False)

    print("  RMSE ratio Ladrillo / BRICK 2.0 (<1: Ladrillo closer to the observations)")
    for wname, _ in WINDOWS:
        w = out[out.window == wname]
        for label, *_ in COMPONENTS:
            l = w[(w.component == label) & (w.arm == "Ladrillo")]
            s = w[(w.component == label) & (w.arm == "BRICK 2.0")]
            if l.empty or s.empty:
                continue
            print(f"    {wname:10s} {label:13s} {l.rmse.iloc[0] / s.rmse.iloc[0]:5.2f}")
    print(f"\nwrote {os.path.relpath(path, REPO) if path.startswith(REPO) else path}")


if __name__ == "__main__":
    main()
