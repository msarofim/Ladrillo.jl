#!/usr/bin/env python3
"""Fig. 4: future sea-level rise to 2300, per component and total, Ladrillo vs BRICK 2.0, on the joint arm.

  python python/figures/plot_future_components.py [--set=vv|ssp] [--arm=joint|fixed]

Writes outputs/figures/future_components_<set>_L27_<arm>.png and its .caption.txt. The paper's Fig. 4 is --set=vv
--arm=joint (the defaults). Reads scripts/project_joint.jl and scripts/brick20_joint.jl products (paths, cells, gates).

PORT of SLR-RFF-BRICK python/plot_future_components.py --set=vv --paper. The drawing code is copied unchanged; only the
I/O differs (figlib.py's header lists how). The research driver also git-logged each scenario's mean-GMST driver to
assert one build; here the forcing files are checked against data/MANIFEST.csv by the package's tests instead.

Both models are on the SAME arm, from the same 2014 splice, the same 1995-2014 re-reference and the same draw-to-config
pairing, so on the joint arm the two bands are the same object and their widths are comparable. They are thinned
differently (Ladrillo 2000 draws before record conditioning, BRICK 2.0 1000; read from the cells files).
"""
import os
import sys
import textwrap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figlib as lf  # noqa: E402

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def _arg(flag, default):
    return next((a[len(flag):] for a in sys.argv[1:] if a.startswith(flag)), default)


SET = _arg("--set=", "vv")
TAG = "L27"                                   # the posterior; part of the output name, as in the research render
ARM = _arg("--arm=", "joint")
if ARM not in ("joint", "fixed"):
    raise SystemExit("--arm must be 'joint' or 'fixed', not %r" % ARM)
DESC = lf.DESC
SCENS = lf.scen_set(SET)
## Above this many scenarios the panel switches from a colour legend to DIRECT LABELS.
DIRECT_LABEL_ABOVE = 3
OUT = lf.fig_path("future_components_%s_%s_%s.png" % (SET, TAG, ARM))

X0, X1 = 2000, 2300
## With three scenarios every band fits; with seven only the two EXTREMES of the set carry one.
BAND_SCENS = ([s[0] for s in SCENS] if len(SCENS) <= 3
              else [SCENS[0][0], SCENS[-1][0]])
ARM_DESC = {"joint": "joint (posterior parameters x 841 FaIR configs)",
            "fixed": "fixed (posterior parameters, mean forcing)"}[ARM]

# --- gates -----------------------------------------------------------------
CHECKS = []
for k, lab, _c, _d in SCENS:
    g = lf.gate_ladrillo(k)
    print("[GATE] %-14s %-7s Ladrillo gates pass (%s, %d configs)" % (lab, k, g["control"], g["configs"]))
    for gate, key, val in g["checks"]:
        CHECKS.append((lab, gate, key, val))
        print("       %s/%s = %+.4f is OVER TOLERANCE (verdict CHECK); stamped on the caption." % (gate, key, val))
CHECK_NOTE = ("" if not CHECKS else
              "  Control over tolerance (drawn anyway, verdict CHECK, not FAIL): "
              + "; ".join("%s %s %+.3f" % (l, k2, v) for l, _g, k2, v in CHECKS)
              + " — a gap between the joint driver's fixed arm and the shipped panel.")


# --- data ------------------------------------------------------------------
def _n_draws(scen, model):
    """The plotted arm's draw count, from the run's cells file (the paths file carries none)."""
    p = lf.require(lf.product_csv(scen, model, "cells"), model, scen)
    d = pd.read_csv(p)
    v = (d[d.arm == ARM] if "arm" in d.columns else d)["n_draws"].unique()
    if len(v) != 1:
        raise SystemExit("[DRAWS] %s carries several draw counts %s" % (os.path.relpath(p, lf.REPO), v))
    return int(v[0])


N_DRAWS = {}
LAD, BRK = {}, {}
for k, lab, _c, _d in SCENS:
    LAD[k] = lf.load_paths(k, "ladrillo", ARM)
    BRK[k] = lf.load_paths(k, "brick20", ARM)
    for _m in ("ladrillo", "brick20"):
        N_DRAWS.setdefault(_m, set()).add(_n_draws(k, _m))
    for name, byc in (("Ladrillo", LAD[k]), ("BRICK 2.0", BRK[k])):
        s, t = lf.check_component_sum(byc, k, name)
        print("[SUM] %-9s %-7s @2300  sum-of-medians %8.2f  vs median-of-total %8.2f  "
              "(%+.2f cm; medians do not add -- reported, not gated)"
              % (name, k, s, t, s - t))

# --- figure ----------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.6))
for ax, comp in zip(axes.ravel(), lf.COMPONENTS):
    ends = []
    for k, lab, col, _d in SCENS:
        d = LAD[k][comp]
        d = d[(d.index >= X0)]
        ax.plot(d.index.values, d.med_cm.values, color=col, lw=1.9, ls="-")
        ends.append((d.med_cm.values[-1], lab, col))
        if k in BAND_SCENS:
            ax.fill_between(d.index.values, d.p05_cm.values, d.p95_cm.values,
                            color=col, alpha=0.13, lw=0)
        e = BRK[k][comp]
        e = e[(e.index >= X0)]
        ax.plot(e.index.values, e.med_cm.values, color=col, lw=1.5, ls="--")
    ax.axvline(2100, color="0.6", lw=0.8, ls=":")
    ax.axhline(0, color="0.85", lw=0.8)
    ax.set_xlim(X0, X1)
    ax.set_title(lf.panel_title(comp), fontsize=10, fontweight="bold", loc="left")
    if len(SCENS) > DIRECT_LABEL_ABOVE:
        lf.direct_label(ax, ends)
    ax.set_ylabel("cm SLE (rel. 1995–2014)", fontsize=8)
    ax.tick_params(labelsize=8)
axes[1, 0].set_xlabel("year")

handles = ([] if len(SCENS) > DIRECT_LABEL_ABOVE
           else [Line2D([], [], color=c, lw=2, label=l) for _k, l, c, _d in SCENS])
handles += [Line2D([], [], color="0.3", ls="-", lw=2, label=DESC["model"]),
            Line2D([], [], color="0.3", ls="--", lw=1.5, label="BRICK 2.0"),
            Line2D([], [], color="0.3", alpha=0.2, lw=8,
                   label="5–95%% band (%s)" % ("all scenarios" if len(SCENS) <= 3
                                               else "extremes only"))]
fig.legend(handles=handles, ncol=min(6, len(handles)), fontsize=8.5, frameon=False,
           loc="upper center", bbox_to_anchor=(0.5, 0.975))
for _m, _v in N_DRAWS.items():
    if len(_v) != 1:
        raise SystemExit("[DRAWS] %s draw count differs across scenarios: %s" % (_m, sorted(_v)))
_NLAD, _NBRK = N_DRAWS["ladrillo"].pop(), N_DRAWS["brick20"].pop()
print("[DRAWS] Ladrillo %d, BRICK 2.0 %d (from the cells files)" % (_NLAD, _NBRK))
_cap = ("%s — %s.  Arm: %s, the same for both models.  %s.  Ladrillo is thinned to %d draws, "
        "BRICK 2.0 to %d.%s"
        % (DESC["model"], DESC["calib"], ARM_DESC, lf.PROJ_BASELINE.capitalize(),
           _NLAD, _NBRK, CHECK_NOTE))
lf.paper_finish(fig, OUT, _cap, rect=[0, 0.10, 1, 0.935])
print("\nwrote %s" % os.path.relpath(OUT, lf.REPO))

# --- console summary: the numbers a caption would quote --------------------
print("\nmedian cm, rel. 1995-2014, %s arm — Ladrillo / BRICK 2.0" % ARM)
for y in (2100, 2150, 2300):
    print("\n  @%d" % y)
    print("    %-14s" % "scenario" + "".join("%22s" % lf.COMP_TITLE[c] for c in lf.COMPONENTS))
    for k, lab, _c, _d in SCENS:
        row = "    %-14s" % lab
        for c in lf.COMPONENTS:
            row += "%11.1f /%9.1f" % (LAD[k][c].med_cm.loc[y], BRK[k][c].med_cm.loc[y])
        print(row)
