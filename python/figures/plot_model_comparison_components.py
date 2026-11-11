#!/usr/bin/env python3
"""Figs. 2 and 3: sea-level rise by component at one horizon, Ladrillo vs BRICK 2.0 vs MAGICC-SLR vs FACTS, on the
seven van Vuuren CMIP7 scenarios.

  python python/figures/plot_model_comparison_components.py [--year=2100|2150|2300|all]

Writes outputs/figures/model_comparison_components_vv_L27_<year>.png and its .caption.txt, one per horizon. The
paper's Fig. 2 is 2100 and Fig. 3 is 2300; 2150 is a sibling. Reads outputs/tables/vv_model_comparison.csv (made by
python/figures/vv_model_comparison.py) and data/comparison/comparator_classes.csv.

PORT of SLR-RFF-BRICK python/plot_model_comparison_components.py --set=vv --year=all --paper. The drawing code is
copied unchanged; only the I/O differs (figlib.py's header lists how). Dropped from the research script: the SSP set
(--set; the package builds no SSP four-source table), --tag, and memo mode (suptitle and in-figure caption). The
research frozen-literature gate is not ported: it compares the SSP set only, so on this set it never ran, and the
caption keeps the sentence the research render stamped in its place (data/comparison/README.md).

THE YEAR IS THE FIGURE, the scenario is the x axis and the component is the panel, so every panel of one figure
shares a horizon and a plausible y-scale. All four bands are drawn: whether a band may be drawn is read from each row's
`band_basis` via figlib.band_is_comparable. FACTS is a stack of modules per component that disagree by up to 8x, so
each module is its own marker; the structured-expert-judgement modules (classified in comparator_classes.csv) get an
open marker, and their 5-95 % whiskers that run past the sampled range are listed with their values, not drawn.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figlib as lf  # noqa: E402

import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def _arg(flag, default):
    return next((a[len(flag):] for a in sys.argv[1:] if a.startswith(flag)), default)


TAG = "L27"                                   # the posterior; part of the output name, as in the research render
DESC = lf.DESC
HORIZONS = [2100, 2150, 2300]
_y = _arg("--year=", "all")
YEARS = HORIZONS if _y == "all" else [int(_y)]
if any(y not in HORIZONS for y in YEARS):
    raise SystemExit("--year must be one of %s or 'all', not %r" % (HORIZONS, _y))

## The scenario set: the research script's --set=vv. Its SET_* tables are kept with the one key the package uses, so
## the drawing code below is the research code unchanged.
SET = "vv"
SET_INFIX = {"vv": "_vv"}
SET_DESC = {"vv": "seven van Vuuren CMIP7 scenarios"}
SET_FIGSIZE = {"vv": (18.0, 8.8)}
SET_XTICK_FS = {"vv": 7.6}
## Past this many off-panel elicited whiskers in one panel, the values move to the caption and the panel keeps a
## one-line pointer. The values are never dropped; only relocated.
CLIP_INPANEL_MAX = 3
## The paper layout keeps the research render's top margin (the figure legend lives there); paper_finish uses only
## the top of this rect.
TOP = 0.925

CMP_CSV = os.path.join(lf.TABLE_DIR, "vv_model_comparison.csv")
CMP_PRODUCER = "python python/figures/vv_model_comparison.py"
CLASSES_CSV = os.path.join(lf.COMPARISON_DIR, "comparator_classes.csv")

## Source order = plotting order = legend order, declared once. `Ladrillo` is the row this
## figure is about and is drawn leftmost in every slot.
SOURCES = ["Ladrillo", "BRICK 2.0", "MAGICC-SLR", "FACTS"]
SRC_LABEL = {"Ladrillo": DESC["model"], "BRICK 2.0": "BRICK 2.0",
             "MAGICC-SLR": "MAGICC-SLR (Nauels 2025)", "FACTS": "FACTS n200 (per module)"}
## Slot centres within one scenario group. FACTS gets the widest slot because it fans its
## modules out inside it; the others are single markers.
SLOT = {"Ladrillo": -0.30, "BRICK 2.0": -0.12, "MAGICC-SLR": 0.06, "FACTS": 0.30}
FACTS_FAN = 0.055        # module-to-module spacing inside the FACTS slot
## OPTION A (Marcus, 2026-09-11): the TOTAL panel draws FACTS as ONE BRACKET, not a fan. The four f-workflows share
## glaciers, TE, LWS and (bar wf4) the GrIS module and differ ONLY in the AIS module, so a fan in the Total panel would
## be the AIS panel's fan drawn a second time. The bracket spans the RANGE OF THE WORKFLOW MEDIANS (thick, capped) and
## the UNION of their 5-95 % (thin); no single "FACTS median" is drawn.
FACTS_BRACKET_PANELS = ("total",)
FACTS_BRACKET_HALF = 0.035   # cap half-width of the bracket
MARK = {"Ladrillo": "s", "BRICK 2.0": "s", "MAGICC-SLR": "D", "FACTS": "o"}
SEJ_MARK = "^"           # open triangle: structured expert judgement, not a model spread

# --- data ------------------------------------------------------------------
if not os.path.exists(CMP_CSV):
    raise SystemExit(
        "no four-source comparison table for %s (set %r) at %s\n  Produce it with:\n    %s"
        % (TAG, SET, os.path.relpath(CMP_CSV, lf.REPO), CMP_PRODUCER))
D = pd.read_csv(CMP_CSV)

missing_src = [s for s in SOURCES if s not in set(D.source)]
if missing_src:
    raise SystemExit("[DATA] %s carries no rows for %s -- this figure is defined by having "
                     "all four sources on one basis, so it must not be drawn with three."
                     % (os.path.relpath(CMP_CSV, lf.REPO), missing_src))

## Ladrillo's module string must BE the posterior being captioned (the table writer stamps it).
_lad_mod = sorted(set(D[D.source == "Ladrillo"].module.astype(str)))
if _lad_mod != [TAG]:
    raise SystemExit("[TAG] %s carries Ladrillo module(s) %s but this figure is captioned "
                     "%s. Refusing to draw." % (os.path.relpath(CMP_CSV, lf.REPO),
                                                _lad_mod, TAG))

# --- gate: the widths this figure DRAWS must be the same object -------------
## THE CAPTION'S COMPARABILITY CLAIM IS DERIVED FROM THE DATA, NEVER TYPED. Every row is asked whether ITS OWN basis
## carries climate uncertainty; the predicate raises on a basis string it does not recognise rather than guessing.
## A cell on a fixed (mean-forcing) arm loses its bar, and the claim weakens, without anyone remembering to.
BASES = {s: sorted(set(D[D.source == s].band_basis.astype(str))) for s in SOURCES}
D["band_ok"] = D.band_basis.map(lf.band_is_comparable)
PARAM_ONLY = sorted({(r.source, b) for r, b in
                     zip(D[~D.band_ok].itertuples(), D[~D.band_ok].band_basis)})
WIDTHS_COMPARABLE = bool(D.band_ok.all())
if WIDTHS_COMPARABLE:
    WIDTH_NOTE = "Bars are 17-83% (thick) and 5-95% (thin), drawn for all four sources."
else:
    WIDTH_NOTE = ("Bars are 17-83%% (thick) and 5-95%% (thin). %d cell(s) are on a "
                  "parameter-only basis and are drawn as medians with no bar: %s."
                  % ((~D.band_ok).sum(),
                     "; ".join("%s = %s" % (a, b) for a, b in PARAM_ONLY)))
for s in SOURCES:
    print("[BASIS] %-11s %s" % (s, " + ".join(BASES[s])))
print("[BASIS] every row carries climate uncertainty: %s"
      % ("YES -- all four bands drawn and comparable" if WIDTHS_COMPARABLE else
         "NO -- %d parameter-only cell(s) drawn WITHOUT a bar: %s"
         % ((~D.band_ok).sum(), PARAM_ONLY)))

# --- the research frozen-literature note -------------------------------------
## The research repository froze its FACTS/MAGICC comparator rows for the SSP set only, so on the van Vuuren set its
## gate stamped this sentence on the caption instead of comparing. Kept verbatim so the caption is the research
## caption; data/comparison/README.md records that this package pins the comparator files by sha256 in
## data/MANIFEST.csv.
LIT_NOTE = ("  the FACTS/MAGICC comparators are frozen for the SSP set only, so "
            "nothing checks that a re-extraction has not moved the %s ones."
            % SET_DESC[SET])
print("[LIT] research frozen arm covers the SSP set only; its caption note is kept on the %s figures "
      "(here the comparator files are pinned by data/MANIFEST.csv)" % SET)

# --- FACTS module classes ---------------------------------------------------
if not os.path.exists(CLASSES_CSV):
    raise SystemExit("[CLASS] no %s -- the structured-expert-judgement modules would then "
                     "be drawn as if they were model spreads. Refusing to draw."
                     % os.path.relpath(CLASSES_CSV, lf.REPO))
_c = pd.read_csv(CLASSES_CSV, comment="#")
CLASS = dict(zip(_c.module.astype(str), _c["class"].astype(str)))
SEJ = sorted(m for m in set(D[D.source == "FACTS"].module.astype(str))
             if CLASS.get(m, "model") == "sej")
print("[CLASS] FACTS modules classified `sej` (open markers; nothing is scored here, the\n        class only changes the marker): %s" % (", ".join(SEJ) or "none"))

# --- console: what each figure can and cannot show --------------------------
print("\n[COVERAGE] sources present per component and horizon")
for y in YEARS:
    for comp in lf.COMPONENTS:
        have = sorted(set(D[(D.year == y) & (D.component == comp)].source))
        gap = [s for s in SOURCES if s not in have]
        print("  @%d  %-20s %-28s %s" % (y, lf.COMP_TITLE[comp], ", ".join(have),
                                         ("MISSING: " + ", ".join(gap)) if gap else ""))


def _sum_check(y):
    """Sum-of-component-medians vs the median of `total`, per source and scenario. REPORTED, NEVER GATED: the sum of
    per-component medians is not the median of the sum unless the components are comonotonic."""
    parts = [c for c in lf.COMPONENTS if c != "total"]
    rows = []
    for src in SOURCES:
        for k, lab, _c2, _d in lf.scen_set(SET):
            s = D[(D.source == src) & (D.scenario == k) & (D.year == y)]
            if src == "FACTS":
                continue          # FACTS has no single module per component to sum
            got = {c: s[s.component == c].med for c in parts + ["total"]}
            if any(len(v) != 1 for v in got.values()):
                continue
            ssum = sum(float(got[c].iloc[0]) for c in parts)
            rows.append((src, lab, ssum, float(got["total"].iloc[0])))
    return rows


# --- figure -----------------------------------------------------------------
SCENS = lf.scen_set(SET)
for YEAR in YEARS:
    OUT = lf.fig_path("model_comparison_components%s_%s_%d.png" % (SET_INFIX[SET], TAG, YEAR))
    absent = sorted({s for s in SOURCES
                     if D[(D.source == s) & (D.year == YEAR)].empty})
    partial = sorted({s for s in SOURCES if s not in absent and
                      set(D[(D.source == s) & (D.year == YEAR)].component)
                      != set(lf.COMPONENTS)})
    CLIP_CAPTION = []          # per-figure; filled by the panels, spent in the caption
    fig, axes = plt.subplots(2, 3, figsize=SET_FIGSIZE[SET])
    for ax, comp in zip(axes.ravel(), lf.COMPONENTS):
        for i, (k, lab, _col, _d) in enumerate(SCENS):
            cell = D[(D.scenario == k) & (D.component == comp) & (D.year == YEAR)]
            for src in SOURCES:
                s = cell[cell.source == src]
                if s.empty:
                    continue
                col = lf.SRC_COLOR[src]
                if src == "FACTS" and comp in FACTS_BRACKET_PANELS:
                    ## OPTION A: one bracket over the workflows (see FACTS_BRACKET_PANELS).
                    mods = sorted(s.module.astype(str))
                    x = i + SLOT[src]
                    m_lo, m_hi = float(s.med.min()), float(s.med.max())
                    ax.plot([x, x], [m_lo, m_hi], color=col, lw=2.6, solid_capstyle="butt",
                            zorder=3)
                    for yv in (m_lo, m_hi):
                        ax.plot([x - FACTS_BRACKET_HALF, x + FACTS_BRACKET_HALF], [yv, yv],
                                color=col, lw=1.6, zorder=3)
                    if s.band_ok.all():
                        ax.plot([x, x], [float(s.p05.min()), float(s.p95.max())], color=col,
                                lw=0.8, alpha=0.7, solid_capstyle="butt", zorder=2)
                    continue
                if src == "FACTS":
                    ## FACTS: one series per MODULE, fanned inside its slot, so the
                    ## DISAGREEMENT between modules is the thing the reader sees.
                    rows = [s[s.module == m].iloc[0] for m in sorted(s.module.astype(str))]
                    x0 = i + SLOT[src] - FACTS_FAN * (len(rows) - 1) / 2
                    xs = [x0 + j * FACTS_FAN for j in range(len(rows))]
                else:
                    rows, xs = [s.iloc[0]], [i + SLOT[src]]
                for r, x in zip(rows, xs):
                    is_sej = CLASS.get(str(r.module), "model") == "sej"
                    mk = SEJ_MARK if is_sej else MARK[src]
                    ms = (5.0 if src == "FACTS" else 6.5)
                    if r.band_ok:
                        ax.plot([x, x], [r.p05, r.p95], color=col,
                                lw=0.8 if src == "FACTS" else 1.0, alpha=0.7,
                                solid_capstyle="butt", zorder=2)
                        ax.errorbar(x, r.med, yerr=[[r.med - r.p17], [r.p83 - r.med]],
                                    fmt=mk, color=col, ms=ms,
                                    capsize=2.0 if src == "FACTS" else 3.5,
                                    lw=1.4 if src == "FACTS" else 2.2,
                                    mfc="none" if is_sej else col, mec=col,
                                    mew=1.4 if is_sej else 0.8, zorder=3)
                    else:
                        ## No bar: a parameter-only width next to three climate-inclusive
                        ## ones would read as a narrow model rather than a narrow arm.
                        ax.plot(x, r.med, mk, ms=ms, color=col,
                                mfc="none" if is_sej else col, mec=col,
                                mew=1.4 if is_sej else 0.8, zorder=3)
        ## Sources absent from THIS panel are named in the panel; a source absent from the WHOLE figure is stamped on
        ## the caption and dropped from the legend instead.
        gone = [s for s in SOURCES
                if s not in absent
                and D[(D.source == s) & (D.component == comp) & (D.year == YEAR)].empty]
        if gone:
            ## A FACTS panel emptied by the climate-driven rule says WHY, not just that.
            _scoped = ("FACTS" in gone and lf.FACTS_CLIMATE_DRIVEN.get(YEAR) is not None)
            _txt = ("FACTS not drawn: no climate-driven module at %d" % YEAR if _scoped and
                    gone == ["FACTS"] else "no " + ", ".join(gone) + " at %d" % YEAR)
            ax.text(0.985, 0.03, _txt,
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=7.2, color="0.35", style="italic")
        ## THE Y-AXIS IS SET BY THE SAMPLED SERIES, AND THE ELICITED ONE IS CLIPPED AND NAMED. The limits come from
        ## the sampling quantiles; any elicited whisker that runs past them is annotated WITH ITS VALUE.
        pan = D[(D.component == comp) & (D.year == YEAR)]
        pan_sej = pan.module.astype(str).map(lambda m: CLASS.get(m, "model") == "sej")
        lo = pd.concat([pan[~pan_sej].p05, pan.med]).min()
        hi = pd.concat([pan[~pan_sej].p95, pan.med]).max()
        pad = 0.08 * max(hi - lo, 1e-6)
        ax.set_ylim(min(lo - pad, 0 if lo >= 0 else lo - pad), hi + pad)
        ## The scenario's DISPLAY label comes from the set table, not from string surgery on the key.
        _lab = {k: l for k, l, _c2, _d in SCENS}
        clipped = ["%s %s %.0f–%.0f" % (_lab.get(r.scenario, r.scenario), r.module,
                                        r.p05, r.p95)
                   for r in pan[pan_sej].itertuples() if r.band_ok and r.p95 > hi + pad]
        if clipped and len(clipped) <= CLIP_INPANEL_MAX:
            ax.text(0.015, 0.975,
                    "elicited 5–95% off panel (cm):\n  " + "\n  ".join(clipped),
                    transform=ax.transAxes, ha="left", va="top", fontsize=6.3,
                    linespacing=1.35, color=lf.SRC_COLOR["FACTS"])
        elif clipped:
            ax.text(0.015, 0.975,
                    "↕ %d elicited 5–95%% ranges off panel — values in caption"
                    % len(clipped),
                    transform=ax.transAxes, ha="left", va="top", fontsize=6.6,
                    color=lf.SRC_COLOR["FACTS"])
            CLIP_CAPTION.append("%s: %s" % (lf.COMP_TITLE[comp], "; ".join(clipped)))
        ax.axhline(0, color="0.85", lw=0.8, zorder=1)
        ax.set_xticks(range(len(SCENS)))
        ax.set_xticklabels([l for _k, l, _c2, _d in SCENS], fontsize=SET_XTICK_FS[SET])
        ax.set_xlim(-0.6, len(SCENS) - 0.4)
        ax.set_title(lf.panel_title(comp), fontsize=10, fontweight="bold", loc="left")
        ax.set_ylabel("cm SLE (rel. 1995–2014)", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(axis="y", alpha=0.25, lw=0.6)

    ## The legend is built from what was actually drawn at this horizon, so a source with
    ## nothing on the figure (e.g. a FACTS module dropped past 2100) gets no entry.
    drawn = [s for s in SOURCES if s not in absent]
    sej_drawn = sorted(set(D[(D.source == "FACTS") & (D.year == YEAR)].module.astype(str))
                       & set(SEJ))
    handles = [Line2D([], [], color=lf.SRC_COLOR[s], marker=MARK[s], ls="none",
                      ms=7, label=SRC_LABEL[s]) for s in drawn]
    if sej_drawn:
        handles += [Line2D([], [], color=lf.SRC_COLOR["FACTS"], marker=SEJ_MARK, ls="none",
                           mfc="none", mew=1.4, ms=7,
                           label="FACTS, structured expert judgement (%s)"
                                 % ", ".join(sej_drawn))]
    _n_wf = D[(D.source == "FACTS") & (D.component == "total") & (D.year == YEAR)].module.nunique()
    if "FACTS" in drawn and _n_wf:
        handles += [Line2D([], [], color=lf.SRC_COLOR["FACTS"], lw=2.6, marker="_", ms=9,
                           mew=1.6, label="FACTS total: range of the %d workflow medians "
                                          "(thick) / union of their 5–95%% (thin)" % _n_wf)]
    handles += [Line2D([], [], color="0.3", lw=2.2,
                       label="17–83% (thick) / 5–95% (thin), all sources"
                             if WIDTHS_COMPARABLE else
                             "17–83% / 5–95% — absent where the arm is parameter-only")]
    fig.legend(handles=handles, ncol=3, fontsize=8.5, frameon=False,
               loc="upper center", bbox_to_anchor=(0.5, 0.972))
    ## CAPTION STYLE (Marcus 2026-09-11b): vintage, baseline, what the bars are, what FACTS
    ## and MAGICC are, and anything absent or unfrozen.
    _keep = lf.FACTS_CLIMATE_DRIVEN.get(YEAR)
    _facts_note = ("the Total panel shows its workflow totals as one bracket over their medians "
                   "(the three emulandice workflows end at 2100; the other four differ only in "
                   "their Antarctic module)" if _keep is None else
                   "only its climate-driven modules are drawn at this horizon (%s); its Greenland "
                   "modules, the other Antarctic modules and every workflow total are not"
                   % ", ".join(sorted(_keep)))
    cap = (
        "%s — %s.  %s.  %s  FACTS n200 per module, rel. baseyear 2005; %s.  "
        "MAGICC-SLR is v7.5.3 + Nauels 2025.  %s%s%s"
        % (DESC["model"], DESC["calib"], lf.PROJ_BASELINE.capitalize(), WIDTH_NOTE, _facts_note,
           ("NOT DRAWN AT %d: %s. " % (YEAR, ", ".join(absent))) if absent else "",
           ("PARTIAL AT %d (some components only): %s. " % (YEAR, ", ".join(partial)))
           if partial else "",
           LIT_NOTE.strip()))
    if CLIP_CAPTION:
        cap += ("  Elicited 5–95%% ranges off panel (cm), listed here rather than over the "
                "data they annotate — %s." % "  ".join(CLIP_CAPTION))
    lf.paper_finish(fig, OUT, cap, rect=[0, 0, 1, TOP])
    plt.close(fig)
    print("\nwrote %s" % os.path.relpath(OUT, lf.REPO))

    print("  [SUM] sum-of-component-medians vs median-of-total @%d "
          "(medians do not add — reported, not gated)" % YEAR)
    for src, lab, ssum, tot in _sum_check(YEAR):
        print("        %-11s %-9s %8.2f vs %8.2f  (%+.2f cm)"
              % (src, lab, ssum, tot, ssum - tot))

# --- console summary: the numbers a caption would quote ---------------------
print("\nmedian cm, rel. 1995–2014.  Ladrillo %s / BRICK 2.0 / MAGICC-SLR / FACTS "
      "[min–max across modules, sej marked *]" % TAG)
for YEAR in YEARS:
    print("\n  @%d" % YEAR)
    for comp in lf.COMPONENTS:
        print("    %s" % lf.COMP_TITLE[comp])
        for k, lab, _c2, _d in SCENS:
            cell = D[(D.scenario == k) & (D.component == comp) & (D.year == YEAR)]

            def one(src):
                s = cell[cell.source == src]
                return "%8.1f" % s.med.iloc[0] if len(s) else "       –"
            f = cell[cell.source == "FACTS"]
            if len(f):
                star = "*" if any(CLASS.get(str(m), "model") == "sej" for m in f.module) else " "
                ftxt = "%7.1f–%-7.1f%s (n=%d)" % (f.med.min(), f.med.max(), star, len(f))
            else:
                ftxt = "            –"
            print("      %-9s %s %s %s  %s"
                  % (lab, one("Ladrillo"), one("BRICK 2.0"), one("MAGICC-SLR"), ftxt))
