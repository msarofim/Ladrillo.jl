#!/usr/bin/env python3
"""Fig. 1: historical sea-level rise 1900-2025 by component and total, Ladrillo vs observations vs BRICK 2.0.

  python python/figures/plot_hindcast_components.py

Writes outputs/figures/hindcast_components_L27.png and its .caption.txt. Reads the posterior-predictive hindcasts
(scripts/posterior_predictive.jl -> outputs/hindcast/postpred_components_timeseries.csv,
scripts/posterior_predictive_brick20.jl -> outputs/hindcast/postpred_brick20_components_timeseries.csv) and the
calibration targets (data/observations/calibration_targets.csv).

PORT of SLR-RFF-BRICK python/plot_hindcast_components.py --paper. The drawing code is copied unchanged; only the I/O
differs (figlib.py's header lists how). The research driver also carried opt-in overlays that the paper figure does not
draw and whose inputs do not ship with this package: MAGICC-SLR's Greenland hindcast (--magicc), the IGCC GMSL line on
the total (--igcc-gmsl), IMBIE 2026 on the ice sheets (--imbie), the >2000 m depth-scope band on thermal expansion
(--deep-band) and a second posterior (--compare). They are dropped here, with their inputs, their startup checks and
their console lines. MAGICC-SLR also entered the baseline gate below; without it the gate's tolerance is set by the
Ladrillo, BRICK 2.0 and observation series alone.

The panels show:
  - every component on one 2x3 grid in figlib.COMPONENTS order, with the land-water panel drawn observation-only
    (neither model emits a land-water hindcast);
  - two Ladrillo bands per component: `_p05/_p95` is posterior-parameter spread, `_pred_p05/_pred_p95` adds the
    calibrated AR(1) + observation error, the band the likelihood asserts.
Both models are integrated from 1850 and saved from 1900 on the same ssp245harm FaIR-mean forcing, so a
Ladrillo-minus-BRICK reading here is a module difference, not a driver difference. The two files spell some components
differently; the name mapping is a declared table below.

BASELINE: everything on this figure is cm rel. 1995-2005, the calibration re-reference, not the 1995-2014 projection
baseline the future figures use. The gate below checks it.

GLACIER OBSERVATION: glaciers are scored against `glaciers_obs_delta_corrected` (the r19-seam adjustment), not the raw
`gsic` target column, so the panel plots that line and re-centres the target's band half-widths on it.
"""
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figlib as lf  # noqa: E402

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

TAG = "L27"                                   # the posterior; part of the output name, as in the research render
DESC = lf.DESC
OUT = lf.fig_path("hindcast_components_%s.png" % TAG)

LAD_CSV = os.path.join(lf.HIND_DIR, "postpred_components_timeseries.csv")
BRK_CSV = os.path.join(lf.HIND_DIR, "postpred_brick20_components_timeseries.csv")
## L27 was fitted to the Frederikse + GRACE build of the targets (ais@1900 = -0.634); the research driver's --targets=
## switch existed to select it over a later IMBIE build. This package ships only that build. Path and md5 are printed.
TGT_CSV = os.path.join(lf.REPO, "data", "observations", "calibration_targets.csv")

BASE0, BASE1 = 1995, 2005          # the CALIBRATION window; see the docstring
X0, X1 = 1900, 2025   # 2025 = last year with AIS/GIS/TE observations (Marcus 2026-10-01)
## X0 is the PLOT span; both arms are integrated from 1850. The caption's span text derives from X0.
## Declared name mapping. Ladrillo/target/BRICK each spell some components differently, and
## a string guess here silently drops a panel.
TGT_COL = {"glaciers": "gsic", "gis": "gis", "ais": "ais", "te": "steric",
           "lws": "lws", "total": "dang"}
BRK_COL = {"glaciers": "gsic", "gis": "gis", "ais": "ais", "te": "te",
           "lws": None, "total": "total"}     # None = no modelled series for that panel
LAD_COL = {c: c for c in lf.COMPONENTS}
LAD_COL["lws"] = None                          # neither model predicts LWS -- panel (e)
                                               # shows the observed series alone
OBS_LINE = {"glaciers": "glaciers_obs_delta_corrected"}   # see the GLACIER OBSERVATION note

C_LAD, C_BRK = lf.SRC_COLOR["Ladrillo"], lf.SRC_COLOR["BRICK 2.0"]
C_OBS = "#333333"
OBS_HATCH = "\\\\\\"                 # back-slash hatch (the research driver's optional >2000 m band used "////")

for f in (LAD_CSV, BRK_CSV, TGT_CSV):
    lf.require(f)
LAD = pd.read_csv(LAD_CSV).set_index("year")
BRK = pd.read_csv(BRK_CSV).set_index("year")
TGT = pd.read_csv(TGT_CSV).set_index("year")
print("[targets] %s  md5 %s  ais@1900 %+.3f" % (TGT_CSV, hashlib.md5(open(TGT_CSV, "rb").read()).hexdigest()[:8],
                                              TGT.loc[1900, "ais"]))


## ---------------------------------------------------------------------------
## BASELINE GATE. Every series on this figure must be on the SAME window, and the two model
## files claim to be already re-referenced to it. That is CHECKED, not trusted: each series'
## own 1995-2005 mean must be ~0. A file silently re-referenced elsewhere would otherwise
## shift a whole panel and look like model bias.
def _base_mean(s):
    w = s.loc[BASE0:BASE1].dropna()
    return float(w.mean()) if len(w) else float("nan")


## The tolerance is derived from what the gate exists to catch. These files are re-referenced
## PER DRAW, so the ensemble MEDIAN over the window need not be exactly zero. The error the gate
## is for is a series baselined to the WRONG WINDOW, so the smallest such displacement present in
## these data (1995-2005 vs the 1995-2014 projection window) sets the scale; the bound is a tenth
## of it, measured here rather than hardcoded.
_alt = []
for _nm, _df, _c in ([("lad", LAD, "%s_p50" % LAD_COL[c]) for c in lf.COMPONENTS if LAD_COL[c]]
                     + [("brk", BRK, "%s_p50" % BRK_COL[c]) for c in lf.COMPONENTS if BRK_COL[c]]
                     + [("obs", TGT, TGT_COL[c]) for c in lf.COMPONENTS]):
    _s = _df[_c]
    _alt.append(abs(_s.loc[BASE0:BASE1].mean() - _s.loc[BASE0:2014].mean()))
BASE_TOL = min(_alt) / 10.0

_off = {}
for comp in lf.COMPONENTS:
    if LAD_COL[comp]:
        _off["Ladrillo/" + comp] = _base_mean(LAD["%s_p50" % LAD_COL[comp]])
    if BRK_COL[comp]:
        _off["BRICK 2.0/" + comp] = _base_mean(BRK["%s_p50" % BRK_COL[comp]])
    _off["obs/" + comp] = _base_mean(TGT[TGT_COL[comp]])
_bad = {k: v for k, v in _off.items() if abs(v) > BASE_TOL}
if _bad:
    raise SystemExit(
        "[BASELINE] these series are NOT zero-mean over %d-%d to %.4f cm, so they are not "
        "on this figure's stated baseline and the panels would be silently offset:\n%s\n"
        "  Do NOT re-reference them here -- fix the producing driver, or the figure's "
        "caption stops being true of its inputs."
        % (BASE0, BASE1, BASE_TOL,
           "\n".join("    %-22s %+.6f cm" % (k, v) for k, v in sorted(_bad.items()))))
print("[BASELINE] all %d model/obs series are zero-mean over %d-%d (max |offset| %.2e cm, "
      "tolerance %.4f cm = 1/10 of the smallest wrong-window displacement in these data)"
      % (len(_off), BASE0, BASE1, max(abs(v) for v in _off.values()), BASE_TOL))

# --- figure ----------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.4))
for ax, comp in zip(axes.ravel(), lf.COMPONENTS):
    t = TGT_COL[comp]
    obs = TGT[t]
    lo, hi = TGT.get("%s_lo" % t), TGT.get("%s_hi" % t)
    if comp == "total":                      # the total target is Dangendorf, +/-1.645 sigma
        lo, hi = TGT["dang_lo"], TGT["dang_hi"]
    if comp in OBS_LINE:                     # glacier seam correction, band re-centred on it
        corr = LAD[OBS_LINE[comp]]
        lo, hi = corr + (lo - obs), corr + (hi - obs)
        obs = corr
    ## HATCHED, NOT FILLED (09-29): a grey obs fill and BRICK 2.0's grey band could not be told
    ## apart where they overlap. Texture, not a second grey, separates them.
    ax.fill_between(obs.index, lo, hi, facecolor="none", edgecolor=C_OBS, hatch=OBS_HATCH,
                    lw=0, alpha=0.40, zorder=1)
    ax.plot(obs.index, obs.values, color=C_OBS, lw=1.6, zorder=4)

    if LAD_COL[comp]:
        c = LAD_COL[comp]
        ax.fill_between(LAD.index, LAD["%s_pred_p05" % c], LAD["%s_pred_p95" % c],
                        color=C_LAD, alpha=0.10, lw=0)
        ax.fill_between(LAD.index, LAD["%s_p05" % c], LAD["%s_p95" % c],
                        color=C_LAD, alpha=0.22, lw=0)
        ax.plot(LAD.index, LAD["%s_p50" % c], color=C_LAD, lw=1.9, zorder=5)
    if BRK_COL[comp]:
        c = BRK_COL[comp]
        ax.fill_between(BRK.index, BRK["%s_p5" % c], BRK["%s_p95" % c],
                        color=C_BRK, alpha=0.16, lw=0)
        ax.plot(BRK.index, BRK["%s_p50" % c], color=C_BRK, lw=1.6, ls="--", zorder=5)
    if not BRK_COL[comp]:
        ## STATED, NOT OMITTED. An empty model panel with no explanation reads as a bug.
        ax.text(0.03, 0.90, "neither model emits an LWS hindcast —\nobservation shown alone",
                transform=ax.transAxes, fontsize=7.6, color="0.35", va="top")

    ax.axhline(0, color="0.85", lw=0.8)
    ax.set_xlim(X0, X1)
    ax.set_title(lf.panel_title(comp), fontsize=10, fontweight="bold", loc="left")
    ax.set_ylabel("cm SLE (rel. 1995–2005)", fontsize=8)
    ax.tick_params(labelsize=8)
axes[1, 0].set_xlabel("year")

## A DRAWN SERIES MUST BE IN THE LEGEND: BRICK 2.0's band and the hatched obs band were both drawn
## but missing from it until 09-29.
handles = [Line2D([], [], color=C_LAD, lw=2, label="%s (median)" % DESC["model"]),
           Patch(facecolor=C_LAD, alpha=0.22, label="%s 5–95%% (parameters)" % "Ladrillo"),
           Patch(facecolor=C_LAD, alpha=0.10, label="%s 5–95%% (predictive, +AR(1)+obs err)" % "Ladrillo"),
           Line2D([], [], color=C_BRK, lw=1.6, ls="--", label="BRICK 2.0 (median)"),
           Patch(facecolor=C_BRK, alpha=0.16, label="BRICK 2.0 5–95%"),
           Line2D([], [], color=C_OBS, lw=1.6, label="observational target"),
           Patch(facecolor="none", edgecolor=C_OBS, hatch=OBS_HATCH, lw=0, alpha=0.40,
                 label="observational target ±1.645σ"),
           ]
## two rows whatever the count: at ncol=2 the legend ran into the top-left panel title
fig.legend(handles=handles, ncol=-(-len(handles) // 2), fontsize=8.5, frameon=False, loc="upper center",
           bbox_to_anchor=(0.5, 0.978))
## CAPTION: what the reader needs to read the panels, plus the provenance labels every output carries.
_cap = (
    ("%s — %s.  Baseline %d–%d.  " % (DESC["model"], DESC["calib"], BASE0, BASE1))
    + "Component observations: Frederikse et al. (2020), 1900–2018, extended by GRACE/GRACE-FO "
    "(AIS, GIS), GlaMBIE 2025 (glaciers), NOAA 0–2000 m thermosteric (TE); total = Dangendorf "
    "2024 extended by NOAA STAR altimetry, and both totals carry the observed land-water storage.  "
    "Both models are run from 1850, plotted from @@X0@@, on the same ssp245harm forcing.  "
    "Observations are the black line, their ±1.645σ range hatched.  ")
_cap = _cap.replace("@@X0@@", str(X0))          # derived from the constant, not retyped
assert "@@" not in _cap, "caption sentinel left unsubstituted"
lf.paper_finish(fig, OUT, _cap, rect=[0, 0.09, 1, 0.915])
print("wrote %s" % os.path.relpath(OUT, lf.REPO))

# --- console summary: 5-year-window comparisons, never single years --------
## A "does it match" check over the historical period is made on at least a 5-year window, so
## interannual variability neither side controls cannot drive the answer. A "1900" label below is
## the 1898-1902 mean.
## MATCHED WINDOWS: every mean is taken over the years where THAT ARM and the obs are BOTH finite,
## and the count is printed. A bare .mean() over the window skips NaN on the obs side only; at the
## ragged modern end (Dangendorf splice ends 2024, GlaMBIE 2023) that compared a 5-year model mean
## with a 3-year obs mean on a rising series and put pure arithmetic into the reported excess.
HALF_WIDTH = 2
EPOCHS = (1900, 1950, 2000, 2024)


def _matched(obs, arm, y):
    """obs and arm means over the years BOTH are finite, plus that count. Never a bare
    .mean() over a window -- that is the bug this function exists to prevent."""
    w = range(y - HALF_WIDTH, y + HALF_WIDTH + 1)
    yy = [t for t in w if t in obs.index and t in arm.index
          and np.isfinite(obs.get(t, np.nan)) and np.isfinite(arm.get(t, np.nan))]
    if not yy:
        return None
    return obs.loc[yy].mean(), arm.loc[yy].mean(), len(yy)


print("\nmodel vs obs, %d-year means centred on each year (cm rel. %d-%d)"
      % (2 * HALF_WIDTH + 1, BASE0, BASE1))
print("  n = years actually used; obs is re-averaged on each arm's own matched window, so the\n"
      "  obs column may differ between arms where their coverage differs. That is the point.")
for y in EPOCHS:
    print("  @%d (%d-%d window)" % (y, y - HALF_WIDTH, y + HALF_WIDTH))
    for comp in lf.COMPONENTS:
        ## THE SAME CORRECTED OBS THE FIGURE PLOTS, so the table cannot contradict the panel.
        _os = LAD[OBS_LINE[comp]] if comp in OBS_LINE else TGT[TGT_COL[comp]]
        row = "    %-22s" % lf.COMP_TITLE[comp]
        for arm_lbl, col_map, df, w in (("Ladrillo", LAD_COL, LAD, 8),
                                        ("BRICK 2.0", BRK_COL, BRK, 8)):
            if not col_map[comp]:
                row += "   %s %*s        " % (arm_lbl, w, "n/a")
                continue
            arm = df["%s_p50" % col_map[comp]]
            r = _matched(_os, arm, y)
            if r is None:
                row += "   %s %*s        " % (arm_lbl, w, "--")
                continue
            o, v, n = r
            row += "   %s obs %*.2f mod %*.2f (%+.2f, n=%d)" % (arm_lbl, w, o, w, v, v - o, n)
        print(row)
