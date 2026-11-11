#!/usr/bin/env python3
"""figlib.py -- shared constants, loaders and gates for the paper's figures (Figs. 1-6).

PORT of SLR-RFF-BRICK python/ladrillo_figs.py, cut down to what the six paper figures use. The constants, the palette,
the panel letters, the direct labelling and the paper-mode layout are copied unchanged, so a figure built here is the
same image as the research render. What differs is the I/O:
  - one posterior (L27), so there is no --tag and no table of calibration vintages;
  - paper mode only: no suptitle and no caption inside the PNG; the caption goes to a sidecar <name>.caption.txt;
  - inputs are this repository's own products (outputs/projections, outputs/hindcast, data/), and the projection gate
    reads the gates file scripts/project_joint.jl writes, which has no CONTROL row (the research driver compared each
    SSP with a shipped panel; this package has none to compare with).
python/figures/check_figures_identity.py compares every PNG with the research render.

This module does not import matplotlib: data and gates only.
"""
import os
import string
import subprocess

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROJ_DIR = os.path.join(REPO, "outputs", "projections")
HIND_DIR = os.path.join(REPO, "outputs", "hindcast")
TABLE_DIR = os.path.join(REPO, "outputs", "tables")
FIG_DIR = os.path.join(REPO, "outputs", "figures")
FORCING_DIR = os.path.join(REPO, "data", "forcing")
COMPARISON_DIR = os.path.join(REPO, "data", "comparison")

# --- the caption's model and calibration (research ladrillo_figs.tag_desc("L27") in paper mode) -----------------
DESC = dict(model="Ladrillo", calib="FaIR 2.2.4 calib 1.6.0 + CMIP7")

# --- baselines -------------------------------------------------------------
PROJ_BASELINE = "cm, rel. 1995-2014"         # every projection product

# --- FACTS: which modules are CLIMATE-DRIVEN at each horizon (Marcus 2026-09-12) ---------
## FittedISMIP GrIS overwrites years >= 2100 with the 2080-2100 rate; bamber19 and deconto21 are lookups keyed on SAT
## integrated over 2000-2099; ar5AIS dynamics are a prescribed function of time; larmip's 200-year response functions
## are zero-padded, so at 2300 every pre-2100 forcing year is dropped. RULE (Marcus): drop any module, and any total
## containing one, whose components are not climate-driven at that horizon (land water excepted).
FACTS_CLIMATE_DRIVEN = {
    2100: None,                                        # None = every module is genuine here
    2150: {"larmip", "ar5glaciers", "tlm", "ssp-lws"},
    2300: {"ar5glaciers", "tlm", "ssp-lws"},
}


def facts_module_ok(module, year):
    """True if this FACTS module/workflow is climate-driven at this horizon (see the table)."""
    keep = FACTS_CLIMATE_DRIVEN.get(int(year))
    return True if keep is None else str(module) in keep


# --- components ------------------------------------------------------------
COMPONENTS = ["glaciers", "gis", "ais", "te", "lws", "total"]
COMP_TITLE = {"glaciers": "Glaciers", "gis": "Greenland ice sheet",
              "ais": "Antarctic ice sheet", "te": "Thermal expansion",
              "lws": "Land-water storage", "total": "Total"}
## Every 2x3 component grid fills its panels in COMPONENTS order, so the letter is a property of the component.
PANEL_LETTER = {c: "(%s)" % string.ascii_lowercase[i] for i, c in enumerate(COMPONENTS)}


def panel_title(comp):
    """'(a) Glaciers': the component's panel letter and its title."""
    return "%s %s" % (PANEL_LETTER[comp], COMP_TITLE[comp])


# --- sources ---------------------------------------------------------------
SRC_COLOR = {"Ladrillo": "#2166ac", "BRICK 2.0": "#7f7f7f",
             "MAGICC-SLR": "#d62728", "FACTS": "#ff9900"}

## Whether a band WIDTH may be drawn is a property of the row's arm (its `band_basis`), not of the source name.
BASIS_CARRIES_CLIMATE = ("joint", "climate + parameter")
BASIS_PARAM_ONLY = ("fixed",)

## Gate verdicts that do not invalidate a run. "ON"/"OFF" record the record-conditioning setting and "onset_year"
## lists a draw the conditioning rejected: records, not checks.
GATE_VERDICTS_OK = ["PASS", "SKIPPED", "measured", "CANCELS", "ON", "OFF", "onset_year"]


def band_is_comparable(basis):
    """True when this row's band carries CLIMATE uncertainty as well as parameter spread."""
    b = str(basis).strip().lower()
    if b.startswith(BASIS_PARAM_ONLY):
        return False
    if b.startswith(BASIS_CARRIES_CLIMATE):
        return True
    raise SystemExit("unrecognised band_basis %r -- refusing to guess whether it carries climate uncertainty."
                     % basis)


# --- scenario sets ---------------------------------------------------------
## (key, label, colour); the van Vuuren table also carries the peak-and-decline flag. The colours were checked for
## colour-vision deficiency in the research repo (python/validate_palette.py); #f69320 is below 3:1 on white, so it is
## always labelled, never the only cue.
SSP_SET = [("ssp126", "SSP1-2.6", "#003466"),
           ("ssp245", "SSP2-4.5", "#f69320"),
           ("ssp585", "SSP5-8.5", "#df0000")]
VV_SET = [("vvVL", "Very Low",      "#00a9cf", True),
          ("vvLN", "Low-to-Neg",    "#1f78b4", True),
          ("vvL",  "Low",           "#003466", False),
          ("vvML", "Medium-to-Low", "#f69320", True),
          ("vvM",  "Medium",        "#c8a000", False),
          ("vvHL", "High-to-Low",   "#df0000", True),
          ("vvH",  "High",          "#7a0002", False)]


def scen_set(name):
    """(key, label, colour, decline) rows for a named set."""
    if name == "ssp":
        return [(k, l, c, False) for k, l, c in SSP_SET]
    if name == "vv":
        return list(VV_SET)
    raise SystemExit("unknown scenario set %r -- expected 'ssp' or 'vv'" % name)


# --- this repository's products ---------------------------------------------
## model 'ladrillo' -> scripts/project_joint.jl; 'brick20' -> scripts/brick20_joint.jl. kind: paths, cells, draws.
PRODUCER = {"ladrillo": ("joint", "scripts/project_joint.jl"), "brick20": ("brick20", "scripts/brick20_joint.jl")}


def product_csv(scen, model, kind):
    prefix, _ = PRODUCER[model]
    return os.path.join(PROJ_DIR, "%s_%s_%s.csv" % (prefix, kind, scen))


def require(path, model=None, scen=None):
    """Stop with the command that makes a missing product."""
    if os.path.exists(path):
        return path
    how = ("julia --project=. %s --scenario=%s" % (PRODUCER[model][1], scen)) if model else "see README.md"
    raise SystemExit("missing %s\n  Produce it with: %s" % (os.path.relpath(path, REPO), how))


def gates_csv(scen):
    return product_csv(scen, "ladrillo", "gates")


def commit_stamp():
    """Short HEAD of this repository, for provenance strings."""
    return subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip() or "UNKNOWN"


# --- gates -----------------------------------------------------------------
## The rows scripts/project_joint.jl writes and this gate requires to be present. A gate that is absent is not a gate
## that passed. (The research driver also wrote CONTROL, OHC-OFFSET, CALIB-MOVE and TRIGGER-PORT rows, which compared
## it with other research products; the package's own tests cover those identities.)
REQUIRED_GATES = ("PAIRING", "SPLICE-IDENTITY", "SUM")


def gate_ladrillo(scen):
    """Read the projection's gates file and refuse a scenario whose run did not pass. Returns the same dict as the
    research gate: control (here: the verdicts of the required gates), configs, and the CHECK rows to stamp."""
    f = require(gates_csv(scen), "ladrillo", scen)
    g = pd.read_csv(f)
    check = g[g.verdict == "CHECK"]
    bad = g[~g.verdict.isin(GATE_VERDICTS_OK + ["CHECK"])]
    if len(bad):
        raise SystemExit("[GATE] %s: %d gate row(s) with an unrecognised verdict:\n%s" % (scen, len(bad), bad))
    missing = [x for x in REQUIRED_GATES if x not in set(g.gate)]
    if missing:
        raise SystemExit("[GATE] %s: gates file has no %s row -- a gate that is absent is not a gate that passed."
                         % (scen, "/".join(missing)))
    _cv = sorted(set(g[g.gate.isin(REQUIRED_GATES)].verdict))
    return dict(control="+".join(_cv),
                configs=int(g.query("gate=='PAIRING' and key=='configs_used'").value.iloc[0]),
                checks=[(r.gate, r.key, float(r.value)) for r in check.itertuples()])


def load_paths(scen, model, arm="joint"):
    """{component: DataFrame indexed by year with med_cm/p05_cm/p95_cm} for one run."""
    f = require(product_csv(scen, model, "paths"), model, scen)
    d = pd.read_csv(f)
    d = d[d.arm == arm]
    if d.empty:
        raise SystemExit("%s %s: arm %r is EMPTY in %s" % (model, scen, arm, os.path.relpath(f, REPO)))
    return {c: g.sort_values("year").set_index("year") for c, g in d.groupby("component")}


def check_component_sum(byc, scen, model, year=2300):
    """Sum of the five component medians vs the median of the total, REPORTED (medians do not add)."""
    parts = [c for c in COMPONENTS if c != "total"]
    missing = [c for c in parts + ["total"] if c not in byc]
    if missing:
        raise SystemExit("[SUM] %s %s: components missing from the trajectory file: %s" % (model, scen, missing))
    s = sum(byc[c].med_cm.loc[year] for c in parts)
    return s, float(byc["total"].med_cm.loc[year])


# --- layout ------------------------------------------------------------------
def direct_label(ax, entries, pad_frac=0.02, fontsize=7.5, min_gap_frac=0.042,
                 collapsed_text="all scenarios"):
    """Label each series at its right-hand end instead of in a legend (Marcus 2026-09-22: the seven van Vuuren
    markers cannot be told apart by colour at the accessibility floors). `entries` is an iterable of (y, text, colour)
    in data coordinates. Labels are pushed apart to at least `min_gap_frac` of the y-range, keeping their order; a set
    that fits inside one gap (land water, identical across scenarios) gets one label."""
    ent = sorted(((float(y), t, c) for y, t, c in entries if y == y), key=lambda e: e[0])
    if not ent:
        return []
    y0, y1 = ax.get_ylim()
    gap = min_gap_frac * (y1 - y0)
    if len(ent) > 1 and (ent[-1][0] - ent[0][0]) < gap:
        ym = sum(e[0] for e in ent) / len(ent)
        x0, x1 = ax.get_xlim()
        ax.set_xlim(x0, x1 + pad_frac * (x1 - x0))
        ax.annotate(collapsed_text, xy=(x1, ym), xytext=(x1 + 0.012 * (x1 - x0), ym),
                    fontsize=fontsize, color="0.35", va="center", ha="left",
                    annotation_clip=False)
        return [ym]
    ys = [e[0] for e in ent]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    for i in range(len(ys) - 2, -1, -1):
        ys[i] = min(ys[i], ys[i + 1] - gap)
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0, x1 + pad_frac * (x1 - x0))
    for (yd, text, col), yl in zip(ent, ys):
        ax.annotate(text, xy=(x1, yd), xytext=(x1 + 0.012 * (x1 - x0), yl),
                    fontsize=fontsize, color=col, va="center", ha="left",
                    annotation_clip=False,
                    arrowprops=None if abs(yl - yd) < 0.25 * gap else
                    dict(arrowstyle="-", color=col, lw=0.5, alpha=0.55,
                         shrinkA=0, shrinkB=0))
    return ys


def fig_path(name):
    """outputs/figures/<name>, the directory created on first use."""
    os.makedirs(FIG_DIR, exist_ok=True)
    return os.path.join(FIG_DIR, name)


def paper_finish(fig, out, caption, rect, dpi=150, **savekw):
    """Lay out and save as the research --paper mode does: no caption band, full canvas, the caller's top margin
    kept (the figure legend lives there), the caption written to <out>.caption.txt."""
    fig.tight_layout(rect=[0, 0.005, 1, rect[3]])
    fig.savefig(out, dpi=dpi, **savekw)
    with open(os.path.splitext(out)[0] + ".caption.txt", "w") as fh:
        fh.write(caption.strip() + "\n")
