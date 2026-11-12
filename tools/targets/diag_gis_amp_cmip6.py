#!/usr/bin/env python3
"""
diag_gis_amp_cmip6.py -- the Greenland amplification shape S(dT) from CMIP6
(what data/constants/gis_amp_shape.csv and its _meta.csv are).

PORT of the research `python/diag_gis_amp_cmip6.py --zone south`
(SLR-RFF-BRICK @ a8a1779, the build that wrote the shipped shape at d59ce12).
load_model, windows, anchor_dt, shape_table and the pooled secant medians are
verbatim. Dropped, because nothing Ladrillo reads depends on them:
  * the per-window table, the per-scenario bootstrap curves and the
    slope-vs-zero bootstrap (their RNG streams never touched the shape: the
    shape is a PCHIP through the pooled MEDIANS, which have no randomness);
  * the summary markdown, the figure (the research figure code raised a
    NameError after the shape was written: `lf` is never imported), the
    `gis_amp_shape_fullcurve` sensitivity arm, and the non-south zones.
Paths come from --sources / --work.

Construction (unchanged): annual CMIP6 tas, global and south-Greenland zone,
anomalies rel 1850-1900 on each model's historical run; rolling 30-yr windows,
dropped if centred before 1950 or below 0.5 K global warming; the SECANT ratio
mean(dT_zone)/mean(dT_glob) per window; median per model within 0.5 K bins,
then median across models (bins with >= 5 models); PCHIP through the bin
medians over 0.75-2.75 K, held flat outside; normalised to 1 at the anchor
dT_eff of the observed fit (diag_gis_amp_anchor.py), inserted as its own grid
node so S(anchor) == 1 exactly under linear interpolation.

The meta file's `commit` column records the commit of the code that built it.
The shipped file says d59ce12 (research); a rebuild here writes the Ladrillo
commit passed as --commit. It is the one cell the identity gate treats as a
provenance label (see check_targets_identity.py).

Reads:  --sources cmip6_gis/tas_series_gis_*.csv
        --work outputs/gis_amp_prior.csv, outputs/gis_amp_anchor.csv
Writes: --work outputs/gis_amp_shape.csv, outputs/gis_amp_shape_meta.csv

Usage: python diag_gis_amp_cmip6.py --sources DIR --work DIR --commit SHORT_SHA
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_t_gis as BTG  # noqa: E402

IN_REL = "cmip6_gis"
OUT_SHAPE_REL = "outputs/gis_amp_shape.csv"
OUT_SHAPE_META_REL = "outputs/gis_amp_shape_meta.csv"
PRIOR_REL = BTG.OUT_PRIORS_REL
ANCHOR_REL = "outputs/gis_amp_anchor.csv"
ANCHOR_SCRIPT = "diag_gis_amp_anchor.py"

# ---- named constants (verbatim) ----------------------------------------------
ZONE = "south"
ZONE_COL = f"tas_gis_{ZONE}"
BASE = BTG.BASE                    # (1850, 1900) anomaly baseline, multi-year rule
SMOOTH_WIN = 30                    # running-window length, years (Antarctic precedent)
YEAR_MIN = 1950                    # pre-1950 windows dropped (small denominator)
DT_MIN = 0.5                       # K; drop windows below this warming level
SCENARIOS = ["ssp126", "ssp245", "ssp585"]
HIST = "historical"
BIN_EDGES = np.arange(0.5, 6.01, 0.5)   # warming-level bins, K
MIN_MODELS_PER_BIN = 5
ESTIMATOR = "secant"
OBS_AMP_FULL_SOUTH_RECEIPT = 1.9221976385152952   # south/full -- the value in use
OBS_AMP_MODERN_SOUTH_RECEIPT = 1.7918384323792236  # south/modern
SHAPE_DT_MIN = 0.75
SHAPE_DT_MAX = 2.75
SHAPE_GRID_STEP = 0.01             # K; the emitted grid the Julia kernel interpolates
SHAPE_GRID_MAX = 8.0               # K; covers ssp585 GMST to 2300

PRIOR_CSV = ANCHOR_CSV = None      # bound in main()


def obs_amp(zone):
    """(full, modern) observed amplification for `zone` from the prior CSV."""
    pr = pd.read_csv(PRIOR_CSV).set_index(["zone", "window"])
    full = float(pr.loc[(zone, "full"), "mean"])
    modern = float(pr.loc[(zone, "modern"), "mean"])
    if zone == "south":
        assert abs(full - OBS_AMP_FULL_SOUTH_RECEIPT) < 1e-9, \
            "south/full prior moved -- gis_amp_prior.csv no longer matches the receipt"
        assert abs(modern - OBS_AMP_MODERN_SOUTH_RECEIPT) < 1e-9, \
            "south/modern prior moved"
    return full, modern


def load_model(path):
    """Return {scenario: DataFrame(year, dt_glob, dt_zone)} anomalies rel BASE."""
    df = pd.read_csv(path)
    if ZONE_COL not in df.columns:
        return None
    df = df.groupby(["scenario", "year"], as_index=False).mean(numeric_only=True)
    hist = df[df.scenario == HIST]
    if hist.empty:
        return None
    seg = hist[(hist.year >= BASE[0]) & (hist.year <= BASE[1])]
    if len(seg) < 20:
        return None
    g0, z0 = seg["tas_global"].mean(), seg[ZONE_COL].mean()

    out = {}
    for sc in SCENARIOS:
        fut = df[df.scenario == sc]
        if fut.empty:
            continue
        joined = (pd.concat([hist, fut])
                    .drop_duplicates(subset="year", keep="last")
                    .sort_values("year"))
        out[sc] = pd.DataFrame({
            "year": joined["year"].values,
            "dt_glob": joined["tas_global"].values - g0,
            "dt_zone": joined[ZONE_COL].values - z0,
        })
    return out


def windows(ts):
    """Rolling SMOOTH_WIN windows -> secant and through-origin slope per window."""
    rows = []
    yr = ts["year"].values
    x, y = ts["dt_glob"].values, ts["dt_zone"].values
    for i in range(len(yr) - SMOOTH_WIN + 1):
        s = slice(i, i + SMOOTH_WIN)
        xc, yc = x[s], y[s]
        xm = xc.mean()
        centre = yr[s].mean()
        if centre < YEAR_MIN or xm < DT_MIN:
            continue
        rows.append({
            "year": centre,
            "dt_glob": xm,
            "secant": yc.mean() / xm,
            "slope": float((xc * yc).sum() / (xc ** 2).sum()),
        })
    return pd.DataFrame(rows)


def anchor_dt():
    """The warming level at which S must equal 1 (diag_gis_amp_anchor.py)."""
    if not os.path.exists(ANCHOR_CSV):
        sys.exit(f"missing {ANCHOR_CSV}; run {ANCHOR_SCRIPT} first")
    a = pd.read_csv(ANCHOR_CSV)
    a = a[(a.zone == ZONE) & (a.window == BTG.AMP_WINDOW_HEADLINE)]
    if a.empty:
        sys.exit(f"{ANCHOR_CSV} has no rows for zone={ZONE} "
                 f"window={BTG.AMP_WINDOW_HEADLINE}")
    return float(a.dt_eff.mean())


def shape_table(pooled, dt_anchor, dt_max=SHAPE_DT_MAX):
    """PCHIP through the pooled binned medians over [SHAPE_DT_MIN, dt_max], held
    flat outside, normalised to 1 at `dt_anchor`."""
    sub = pooled[(pooled.dt_bin >= SHAPE_DT_MIN) & (pooled.dt_bin <= dt_max)]
    pch = PchipInterpolator(sub["dt_bin"].values, sub["median"].values,
                            extrapolate=False)
    clip = lambda d: np.clip(d, SHAPE_DT_MIN, dt_max)  # noqa: E731
    r_anchor = float(pch(clip(dt_anchor)))
    dt = np.round(np.arange(0.0, SHAPE_GRID_MAX + SHAPE_GRID_STEP / 2,
                            SHAPE_GRID_STEP), 4)
    dt = np.unique(np.concatenate([dt, [dt_anchor]]))
    r = pch(clip(dt))
    return pd.DataFrame({"dt": dt, "R_secant": r, "S": r / r_anchor}), r_anchor


def main():
    global PRIOR_CSV, ANCHOR_CSV
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sources", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--commit", required=True, help="short commit of the code that builds this (meta column)")
    a = ap.parse_args()
    PRIOR_CSV = os.path.join(a.work, PRIOR_REL)
    ANCHOR_CSV = os.path.join(a.work, ANCHOR_REL)
    OUT_SHAPE = os.path.join(a.work, OUT_SHAPE_REL)
    OUT_SHAPE_META = os.path.join(a.work, OUT_SHAPE_META_REL)
    OBS_AMP_FULL, _ = obs_amp(ZONE)

    files = sorted(glob.glob(os.path.join(a.sources, IN_REL, "tas_series_gis_*.csv")))
    if not files:
        sys.exit(f"no inputs in {os.path.join(a.sources, IN_REL)}")

    recs, skipped = [], []
    for f in files:
        model = os.path.basename(f)[len("tas_series_gis_"):-len(".csv")]
        try:
            per_sc = load_model(f)
        except Exception as e:
            skipped.append(f"{model}: {type(e).__name__}: {e}")
            continue
        if not per_sc:
            skipped.append(f"{model}: no usable scenario/baseline")
            continue
        for sc, ts in per_sc.items():
            w = windows(ts)
            if w.empty:
                continue
            w["model"], w["scenario"] = model, sc
            recs.append(w)
    if not recs:
        sys.exit("no usable windows")
    allw = pd.concat(recs, ignore_index=True)
    n_models = allw.model.nunique()

    # ---- pooled (all-scenario) secant curve: the deliverable shape ------------
    sub = allw.copy()
    sub["bin"] = pd.cut(sub["dt_glob"], BIN_EDGES,
                        labels=(BIN_EDGES[:-1] + BIN_EDGES[1:]) / 2)
    pooled = (sub.groupby(["bin", "model"], observed=True)[ESTIMATOR].median()
                 .groupby("bin", observed=True)
                 .agg(["median", "count"]).reset_index())
    pooled = pooled[pooled["count"] >= MIN_MODELS_PER_BIN]
    pooled["dt_bin"] = pooled["bin"].astype(float)

    dt_anchor = anchor_dt()
    grid, r_anchor = shape_table(pooled, dt_anchor)
    # full precision, no float_format: the anchor node must round-trip exactly
    os.makedirs(os.path.dirname(OUT_SHAPE), exist_ok=True)
    grid.to_csv(OUT_SHAPE, index=False)
    pd.DataFrame([dict(anchor_dt=dt_anchor, r_anchor=r_anchor,
                       obs_amp_full=OBS_AMP_FULL, estimator=ESTIMATOR,
                       dt_min=SHAPE_DT_MIN, dt_max=SHAPE_DT_MAX,
                       grid_step=SHAPE_GRID_STEP, grid_max=SHAPE_GRID_MAX,
                       n_models=n_models, zone=ZONE, commit=a.commit)]
                ).to_csv(OUT_SHAPE_META, index=False)
    print(f"diag_gis_amp_cmip6 | {n_models} models, {len(allw)} windows"
          f"{'; skipped ' + '; '.join(skipped) if skipped else ''} | anchor {dt_anchor:.4f} K, "
          f"R_secant(anchor) {r_anchor:.4f} | wrote {OUT_SHAPE_REL}, {OUT_SHAPE_META_REL}")


if __name__ == "__main__":
    main()
