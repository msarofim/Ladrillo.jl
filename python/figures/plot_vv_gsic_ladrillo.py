#!/usr/bin/env python3
"""Fig. 5: glacier melt to 2300 under the seven van Vuuren CMIP7 scenarios, Ladrillo's 3-reservoir nu glaciers.

  python python/figures/plot_vv_gsic_ladrillo.py [--arm=joint|fixed]

Writes outputs/figures/vv_gsic_ladrillo_L27_2300_<arm>.png (the fixed render drops the arm suffix, as in the research
repo) and its .caption.txt. The paper's Fig. 5 is --arm=joint, the default. Panels: (a) the scenario GMST forcing,
(b) cumulative glacier melt with the 5-95% band on the Very Low scenario, (c) the melt rate on the four
peak-and-decline scenarios.

PORT of SLR-RFF-BRICK python/plot_vv_gsic_wr_vs_ladrillo.py --ladrillo-only --paper. The drawing code is copied
unchanged; only the I/O differs:
  - --ladrillo-only and --paper are resolved: the Wigley-Raper (BRICK 2.0) panel, its lines and its legend entries are
    removed (with them the melt-rate panel's legend, which is empty under --ladrillo-only), and there is no title and
    no in-figure caption.
  - Panel (a) reads data/forcing/fair_mean_gmst_<marker>.csv. The research script drew panel (a) from the `gmst` column
    of outputs/vv_gsic_2300.csv (the Wigley-Raper arm's product) and gated that column to equal these driver files
    within 1e-6 K; over the drawn years 2000-2300 the two are equal to the bit. The Wigley-Raper file is not read.
  - Glacier paths and gates come from scripts/project_joint.jl (outputs/projections/joint_{paths,gates}_<marker>.csv).
    figlib.gate_ladrillo requires the PAIRING, SPLICE-IDENTITY and SUM rows (there is no CONTROL row to require); as in
    the research script, a CHECK verdict also stops the figure.
  - The research script git-logged each scenario's mean-GMST driver to assert one build and put that commit in the
    caption. Here the forcing files are checked against data/MANIFEST.csv by the package's tests instead, and the
    caption says so.

Under the joint arm each posterior draw runs on its own FaIR configuration, so the band carries parameter plus climate
uncertainty; panel (a) draws the ensemble-mean GMST and says so. Ladrillo maps that GMST to each glacier block's
temperature as amp_b * GMST, spliced to the observed record, so the panel shows the forcing, not the temperature each
reservoir sees.
"""
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


def _arg(flag, default):
    return next((a[len(flag):] for a in sys.argv[1:] if a.startswith(flag)), default)


LADRILLO_TAG = "L27"                          # the posterior; part of the output name, as in the research render
ARM = _arg("--arm=", "joint")
if ARM not in ("fixed", "joint"):
    raise SystemExit("--arm must be fixed or joint, not %r" % ARM)
OUTPNG = lf.fig_path("vv_gsic_ladrillo_%s_2300%s.png" % (LADRILLO_TAG, "" if ARM == "fixed" else "_" + ARM))
GMST_CSV = os.path.join(lf.FORCING_DIR, "fair_mean_gmst_{m}.csv")

## lf.VV_SET rows are (key, label, colour, decline); MARKERS is (label, key, colour, decline).
MARKERS = [(lab, k, c, d) for k, lab, c, d in lf.VV_SET]
LABELS = [m[0] for m in MARKERS]
KEY = {m[0]: m[1] for m in MARKERS}
COL = {m[0]: m[2] for m in MARKERS}
DECLINE = [m[0] for m in MARKERS if m[3]]
SPREAD_LO, SPREAD_HI = "Very Low", "High"
X0, X1 = 2000, 2300
## The paper names the model by its version, not the calibration tag; the tag stays in the filename and the caption.
LAD_NAME = "Ladrillo (3-reservoir ν)"

# --- gates and data ----------------------------------------------------------
## A CHECK verdict stops the figure: the research script accepted only lf.GATE_VERDICTS_OK.
LAD, GMST = {}, {}
_nc = None
for lab in LABELS:
    m = KEY[lab]
    g = lf.gate_ladrillo(m)
    if g["checks"]:
        raise SystemExit("[GATE] %s: %d gate row(s) with verdict CHECK:\n%s" % (m, len(g["checks"]), g["checks"]))
    _nc = g["configs"] if _nc is None else _nc
    byc = lf.load_paths(m, "ladrillo", ARM)
    if "glaciers" not in byc:
        raise SystemExit("[GATE] %s: no `glaciers`/`%s` rows in %s"
                         % (m, ARM, os.path.relpath(lf.product_csv(m, "ladrillo", "paths"), lf.REPO)))
    LAD[lab] = byc["glaciers"]
    gf = lf.require(GMST_CSV.format(m=m))
    t = pd.read_csv(gf).rename(columns={"gmst_C": "gmst"}).sort_values("year")
    t = t[(t.year >= X0) & (t.year <= X1)][["year", "gmst"]]
    ## The forcing and the projection must cover the same drawn years (the driver files run one year past 2300).
    _ly = LAD[lab].index.values
    _ly = _ly[_ly >= X0]
    if not np.array_equal(t.year.values, _ly):
        raise SystemExit("[YEARS] %s: from %d, %s covers %d-%d and the glacier paths %d-%d"
                         % (m, X0, os.path.relpath(gf, lf.REPO), t.year.min(), t.year.max(), _ly.min(), _ly.max()))
    GMST[lab] = t
print("[GATE] all %d Ladrillo %s arms pass (arm=%s, %d configs)" % (len(LABELS), LADRILLO_TAG, ARM, _nc))


def gmst(lab):
    return GMST[lab]


def lad(lab):
    d = LAD[lab]
    return d[d.index >= X0]


fig, _axs = plt.subplots(3, 1, figsize=(8.6, 9.4), sharex=True,
                         gridspec_kw=dict(height_ratios=[0.85, 1.2, 0.95], hspace=0.13))
ax = {0: _axs[0], 1: None, 2: _axs[1], 3: _axs[2]}   # same panel keys as the research script, WR slot empty

# ---- (a) GMST forcing ----
## Direct labels, not a colour legend (Marcus's ruling 2026-09-22): see figlib.direct_label.
_ends0 = []
for s in LABELS:
    d = gmst(s)
    ax[0].plot(d.year.values, d.gmst.values, color=COL[s], lw=1.8)
    _ends0.append((d.gmst.values[-1], s, COL[s]))
ax[0].set_ylabel("GMST (°C rel. PI)")
ax[0].text(0.012, 0.93, "(a)  GMST forcing" + (" (FaIR ensemble mean)" if ARM == "joint" else ""),
           transform=ax[0].transAxes, fontsize=9.5, fontweight="bold", va="top")
lf.direct_label(ax[0], _ends0)
_pk = {s: (gmst(s).set_index("year").gmst.loc[2015:2300].idxmax(),
           gmst(s).set_index("year").gmst.loc[2015:2300].max()) for s in DECLINE}
ax[0].annotate("%d peak-and-decline pathways\n(peaks %.2f–%.2f °C, %d–%d)"
               % (len(DECLINE), min(v for _, v in _pk.values()),
                  max(v for _, v in _pk.values()),
                  min(y for y, _ in _pk.values()), max(y for y, _ in _pk.values())),
               xy=(2235, 2.55), fontsize=7.5, color="0.3", ha="center")

_lad_max = max(lad(s).p95_cm.max() for s in LABELS)
ymax = _lad_max * 1.02

# ---- (b) Ladrillo ----
_ends2 = []
for s in LABELS:
    d = lad(s)
    ax[2].plot(d.index.values, d.med_cm.values, color=COL[s], lw=1.9)
    _ends2.append((d.med_cm.values[-1], s, COL[s]))
    if s == SPREAD_LO:
        ax[2].fill_between(d.index.values, d.p05_cm.values, d.p95_cm.values,
                           color=COL[s], alpha=0.15, lw=0)
ax[2].text(0.012, 0.93, "(b)" + "  %s — the 3 reservoirs equilibrate, so declining T slows the melt"
           % LAD_NAME, transform=ax[2].transAxes, fontsize=9.5, fontweight="bold", va="top")
sp_lad = lad(SPREAD_HI).med_cm.loc[2300] - lad(SPREAD_LO).med_cm.loc[2300]
ax[2].legend(handles=[Line2D([], [], color="none",
                             label=f"{SPREAD_LO}→{SPREAD_HI} spread @2300:")]
                     + [Line2D([], [], color="none",
                               label=f"   {LAD_NAME}  {sp_lad:.1f} cm")],
             fontsize=8, loc="lower right", handlelength=0, handletextpad=0,
             frameon=True, framealpha=0.88, edgecolor="none", facecolor="white")

for a in [x for x in (ax[1], ax[2]) if x is not None]:
    a.set_ylim(0, ymax)
    a.axvline(2100, color="0.6", lw=0.8, ls=":")
    a.set_ylabel("cumulative glacier\nmelt (cm SLE, rel 1995–2014)")
## After set_ylim: direct_label derives its minimum label gap from the y-range.
lf.direct_label(ax[2], _ends2)

# ---- (c) melt rate on the four decline pathways ----
## Built from the DECLINE flag in MARKERS, so it cannot fall out of step with panel (a).
_ends3 = []
for s in DECLINE:
    e = lad(s)
    _r = np.gradient(e.med_cm.values, e.index.values) * 100
    ax[3].plot(e.index.values, _r, color=COL[s], lw=1.8, ls="-")
    _ends3.append((_r[-1], s, COL[s]))
ax[3].axhline(0, color="0.6", lw=0.8)
ax[3].axvline(2100, color="0.6", lw=0.8, ls=":")
ax[3].set_ylabel("melt rate\n(cm / century)")
ax[3].set_xlabel("year")
ax[3].set_xlim(X0, X1)
lf.direct_label(ax[3], _ends3)
ax[3].text(0.012, 0.93, ("(c)  melt rate on the %d peak-and-decline pathways — falls toward zero "
                         "as the reservoirs equilibrate") % len(DECLINE),
           transform=ax[3].transAxes, fontsize=9.5, fontweight="bold", va="top")
## The research script also drew ax[3].legend(handles=[]) here (its handles are the WR/Ladrillo line styles, empty
## under --ladrillo-only). It draws nothing, and the PNG is byte-identical without it, so it is not ported.

## Caption scope: whose posterior, which arm, which forcing.
_SPREAD = ("each posterior draw on its own FaIR configuration, so the band is parameter plus climate "
           "uncertainty; panel (a) is the ensemble-mean GMST" if ARM == "joint" else
           "posterior-parameter spread on mean forcing")
_CAPTION = (("Ladrillo %s 3-reservoir ν glaciers (%s chains, `%s` arm), " + _SPREAD + ".\n")
            % (LADRILLO_TAG, LADRILLO_TAG, ARM)
            + "FaIR 2.2.4 (calib 1.6.0) van Vuuren scenario GMST; one build, one calibration "
              "throughout (the forcing files are checked against data/MANIFEST.csv by the package's tests), "
              "each scenario on its own CMIP7 land-use, irrigation and volcanic/solar forcing.")
with open(os.path.splitext(OUTPNG)[0] + ".caption.txt", "w") as _fh:
    _fh.write(_CAPTION.replace("\n", " ").strip() + "\n")
fig.savefig(OUTPNG, dpi=150, bbox_inches="tight")
print("wrote " + os.path.relpath(OUTPNG, lf.REPO))

print("\n%-16s %8s %8s %10s %12s" % ("scenario", "peak K", "@yr", "Lad@2300", "Lad rate@2300"))
for s in LABELS:
    g = gmst(s).set_index("year").gmst.loc[2015:2300]
    e = lad(s)
    print("%-16s %8.2f %8d %10.2f %12.2f"
          % (s, g.max(), g.idxmax(), e.med_cm.loc[2300],
             np.gradient(e.med_cm.values, e.index.values)[-1] * 100))
