#!/usr/bin/env python3
"""build_shared_climate_nc.py — FACTS climate NetCDFs on the SHARED climate-driver convention.

  python3 build_shared_climate_nc.py --scen=ssp245 [--nsamp=200]
  python3 build_shared_climate_nc.py --all

WHAT "SHARED" MEANS, AND WHY THIS SUPERSEDES build_fair_climate_nc.py.
That script was a PoC: hardcoded to ssp245, to calib 1.4.5, to one harmonized CSV, and it
fed FACTS the RAW per-config path. Ladrillo and BRICK 2.0 are not driven by the raw path --
`julia/scope_slr_fair_uncertainty.jl` applies the SPLICED convention (its `convention()`,
~:226). Handing FACTS a different path from the other two models would make "one shared
climate driver" false at the moment of writing it down.

So this reproduces the Julia convention EXACTLY, in Python:
    y <= 2014 : the shipped MEAN driver, untouched  (so the calibrated hindcast is identical)
    y >  2014 : mean(REF) + (config path - mean(config over REF)),  REF = 1995-2014
Uncertainty enters the FUTURE only. [SPLICE-MATCH] below proves the reproduction rather
than asserting it, against the same committed cube and the same 2e-6 degC tolerance the
Julia driver gates on.

THE CONVENTION THIS ESTABLISHES, to be stated on every product built from these files:
  Ladrillo / BRICK 2.0 / FACTS       -> ONE climate: FaIR 2.2.4 calib 1.6.0 + CMIP7,
                                        the same 841-config cubes, the same 2014 splice,
                                        the same 1995-2014 reference, on BOTH scenario sets.
  MAGICC-SLR                         -> its OWN climate, from van Vuuren / SSP EMISSIONS.
                                        That is not a gap in the convention, it is the
                                        point of keeping one independent arm.

⚠ THE 2300->2500 TAIL IS A CONSTANT HOLD AND MUST NEVER REACH A REPORTED YEAR. FACTS wants a
1750-2500 axis; our cubes stop at 2300, which is the standing end year for all our work. The
hold is a padding artefact, so [HORIZON] refuses to build if FACTS could report past 2300.
"""
import os
import sys

import numpy as np
import netCDF4 as nc4
import pandas as pd

OBS = os.path.expanduser(
    "~/Documents/2026/CodeProjects/SLR-RFF-BRICK/data/observations")
FACTS_REPO = os.path.dirname(os.path.abspath(__file__))

## The two scenario sets, declared.
##
## ⚠ THE SCENARIO STRING IS A GROUP KEY, NOT A PHYSICS INPUT -- corrected 2026-08-31 after
## Marcus asked what the "SSP host" involved. An earlier draft here keyed each van Vuuren
## scenario to a host SSP (vvH -> ssp370, ...). That was wrong-headed: the driving path is
##      van Vuuren emissions -> FaIR 2.2.4 calib 1.6.0 (841 cfgs) -> GMST + OHC -> FACTS
## and NOTHING on it consults an SSP. The string only names the group a module reads its
## climate from (`/<scen>/surface_temperature`), and we own both sides of that lookup, so
## the scenario's own key is the honest name. Every module that matters here is driven by the
## injected CLIMATE -- including bamber19, which in injected mode skips its scenario_map
## branch entirely and picks its high/low core from the integrated SAT over 2000-2099 of the
## climate we hand it (`pickScenario`).
##
## THE THREE PLACES A SCENARIO STRING STILL REACHES A NUMBER, all checked in the sources:
##   * ipccar5 icesheets has `if scenario in ['rcp85','ssp585']`, which widens the GREENLAND
##     dynamics range -- but our config runs that module through `pipeline.AIS.global.yml`,
##     which emits only AIS/EAIS/WAIS. INERT for this comparison.
##   * emulandice reads the `Scenario` attribute to select per-SSP-TRAINED GP emulators.
##     Genuinely label-dependent, and 2100-capped. emulandice runs on the ten base keys
##     (three SSP controls + seven van Vuuren scenarios) at 2100 (2026-09-14).
##   * ssp/landwaterstorage takes scenario=ssp1..ssp5 to pick a POPULATION pathway
##     (reservoir construction + groundwater depletion). This one is real, and it is not
##     climate at all -- it cannot be derived from temperature. But the van Vuuren scenarios
##     genuinely CARRY an SSP socioeconomic basis; it is in their own filenames. So `lws_ssp`
##     below is the scenario's OWN socioeconomics, not a climate proxy.
SCENARIOS = [
    # key      scen key   lws_ssp  label
    ("ssp126", "ssp126",  "ssp1", "SSP1-2.6"),
    ("ssp245", "ssp245",  "ssp2", "SSP2-4.5"),
    ("ssp585", "ssp585",  "ssp5", "SSP5-8.5"),
    ("vvVL",   "vvVL",    "ssp1", "van Vuuren Very Low scenario (very_low_ssp1_marker)"),
    ("vvLN",   "vvLN",    "ssp2", "van Vuuren Low-to-Negative scenario (low_to_negative_ssp2_marker)"),
    ("vvL",    "vvL",     "ssp2", "van Vuuren Low scenario (low_ssp2_marker)"),
    ("vvML",   "vvML",    "ssp2", "van Vuuren Medium-to-Low scenario (medium_to_low_ssp2_marker)"),
    ("vvM",    "vvM",     "ssp2", "van Vuuren Medium scenario (medium_ssp2_marker)"),
    ("vvHL",   "vvHL",    "ssp5", "van Vuuren High-to-Low scenario (high_to_low_ssp5_marker)"),
    ("vvH",    "vvH",     "ssp3", "van Vuuren High scenario (high_ssp3_marker)"),
    ## THE MATCHED-dT OVERSHOOT PAIR (2026-09-02). Our native ssp534over crosses BELOW
    ## ssp126 in 2127 and ends 0.126 K COOLER, so the native pair cannot answer SLEIP's
    ## overshoot question at all. `ssp534overMATCH` = ERF_126 + max(ERF_534-ERF_126, 0),
    ## built by FaIRtoFrEDI scripts/build_fair_cube_matched_dt.py; GSAT relaxes back FROM
    ## ABOVE (+0.042 K @2150). ⚠ IDEALISED -- never label it SSP5-3.4-OS.
    ## ⭐ BOTH ARMS TAKE lws_ssp="ssp1", DELIBERATELY. landwaterstorage picks a POPULATION
    ## pathway, which is not climate and does not belong in an overshoot penalty. Giving
    ## the overshoot arm ssp5 would put a reservoir/groundwater difference INTO the
    ## penalty. Holding it fixed makes LWS cancel, exactly as it does in Ladrillo and
    ## BRICK 2.0, where the measured LWS penalty is 0.000 cm at every horizon.
    ("ssp126nomarker",   "ssp126nomarker",   "ssp1", "SSP1-2.6 marker-free (penalty reference)"),
    ("ssp534overMATCH",  "ssp534overMATCH",  "ssp1", "matched-dT overshoot arm (IDEALISED, not SSP5-3.4-OS)"),
    ## THE SAME PAIR AT 2300 (2026-09-02). SLEIP's overshoot headline is a 2300 number and our
    ## 2150 arms could not reach it. Separate KEYS, not a change to the ones above: the existing
    ## arms keep pyear_end 2150 so the shipped FACTS comparison column stays like-for-like, and
    ## these two are additive. ⚠ 2300 is UNTESTED for this module set — that is why it is a new
    ## arm rather than an edit to an old one.
    ("ssp126nomarker2300",  "ssp126nomarker2300",  "ssp1", "SSP1-2.6 marker-free, to 2300"),
    ("ssp534overMATCH2300", "ssp534overMATCH2300", "ssp1", "matched-dT overshoot arm, to 2300 (IDEALISED)"),
    ## THE SEVEN VAN VUUREN SCENARIOS AT 2300 (2026-09-03), for the deliverable's Fig 4 (the
    ## 2300 van Vuuren component comparison), which had no FACTS comparator at 2300 at all --
    ## only the matched-dT pair above went past 2150. Same pattern: separate keys, existing
    ## 2150 arms untouched so the shipped comparison stays like-for-like.
    ("vvVL2300", "vvVL2300", "ssp1", "van Vuuren Very Low scenario (very_low_ssp1_marker), to 2300"),
    ("vvLN2300", "vvLN2300", "ssp2", "van Vuuren Low-to-Negative scenario (low_to_negative_ssp2_marker), to 2300"),
    ("vvL2300",  "vvL2300",  "ssp2", "van Vuuren Low scenario (low_ssp2_marker), to 2300"),
    ("vvML2300", "vvML2300", "ssp2", "van Vuuren Medium-to-Low scenario (medium_to_low_ssp2_marker), to 2300"),
    ("vvM2300",  "vvM2300",  "ssp2", "van Vuuren Medium scenario (medium_ssp2_marker), to 2300"),
    ("vvHL2300", "vvHL2300", "ssp5", "van Vuuren High-to-Low scenario (high_to_low_ssp5_marker), to 2300"),
    ("vvH2300",  "vvH2300",  "ssp3", "van Vuuren High scenario (high_ssp3_marker), to 2300"),
    ## THE THREE SSP CONTROLS AT 2300 (2026-09-11, Marcus's 9/11 comment [12]): the SSP
    ## comparison figure at 2300 (deliverable FIG 8) had no FACTS bars because only the van
    ## Vuuren set was ever extended. Same pattern again: separate keys, the 2150 controls
    ## untouched. emulandice runs on the ten base keys at 2100 (2026-09-14), not on these
    ## *2300 twins: it is 2100-capped (nothing to add at 2300) and selects per-SSP-TRAINED
    ## emulators from the Scenario label, which here is the key itself.
    ## [EXTEND] in extract_facts_shared_components.py requires each *2300 twin to reproduce
    ## its base run BIT-FOR-BIT on the overlapping years before its tail is merged.
    ("ssp1262300", "ssp1262300", "ssp1", "SSP1-2.6, to 2300"),
    ("ssp2452300", "ssp2452300", "ssp2", "SSP2-4.5, to 2300"),
    ("ssp5852300", "ssp5852300", "ssp5", "SSP5-8.5, to 2300"),
]
SCEN_BY_KEY = {k: (f, lws, l) for k, f, lws, l in SCENARIOS}

## ⭐ THE EXPERIMENT KEY AND THE CUBE FILENAME ARE NOT THE SAME STRING, and they cannot be.
## RADICAL EnTK refuses a pipeline name containing anything but letters, digits and '.'
## (entk/pipeline.py:169), so an UNDERSCORE in the key kills the run one second in with a
## message about `pipeline.0000`, nowhere near the cause. Our cube convention is
## `fair_cube_gmst_<ssp>_nomarker_raw.csv` -- underscored by construction -- so every
## marker-free scenario would hit this. The key stays EnTK-safe; CUBE_KEY maps it back to
## the file on disk. Measured 2026-09-02: `ssp534overMATCH` ran fine, `ssp126_nomarker`
## failed in 1 s, which is exactly the difference.
CUBE_KEY = {"ssp126nomarker": "ssp126_nomarker",
            "ssp126nomarker2300": "ssp126_nomarker",
            "ssp534overMATCH2300": "ssp534overMATCH",
            "vvVL2300": "vvVL", "vvLN2300": "vvLN", "vvL2300": "vvL", "vvML2300": "vvML",
            "vvM2300": "vvM", "vvHL2300": "vvHL", "vvH2300": "vvH",
            "ssp1262300": "ssp126", "ssp2452300": "ssp245", "ssp5852300": "ssp585"}

## Per-key reporting horizon. Default is the 2150 convention every prior arm used; the `*2300`
## keys are the only ones that go further. The value TRAVELS IN THE NetCDF (attribute `PyearEnd`)
## so build_shared_configs.py reads it from the file rather than re-typing it — the same rule the
## Scenario and LWSScenario fields already follow, so a config cannot disagree with its climate.
HORIZON_BY_KEY = {"ssp126nomarker2300": 2300, "ssp534overMATCH2300": 2300,
                  "vvVL2300": 2300, "vvLN2300": 2300, "vvL2300": 2300, "vvML2300": 2300,
                  "vvM2300": 2300, "vvHL2300": 2300, "vvH2300": 2300,
                  "ssp1262300": 2300, "ssp2452300": 2300, "ssp5852300": 2300}
def horizon(key):
    return HORIZON_BY_KEY.get(key, FACTS_REPORT_MAX)

## The guard that turns that cryptic failure into a refusal AT BUILD TIME.
_BAD = [k for k, _f, _w, _l in SCENARIOS if not __import__("re").fullmatch(r"[A-Za-z0-9.]+", k)]
if _BAD:
    raise SystemExit("[KEY] EnTK forbids anything but letters/digits/'.' in a pipeline name; "
                     "these keys would fail one second into the run: %s" % _BAD)

SPLICE_YEAR = 2014
REF = (1995, 2014)                 # LADRILLO_REF, ladrillo_projection.jl:255
SPLICE_MATCH_TOL = 2e-6            # degC; scope_slr_fair_uncertainty.jl:114
FACTS_YEARS = np.arange(1750, 2501)
## ⚠ THIS IS OUR CONVENTION, NOT A FACTS LIMIT. Corrected 2026-09-02 after Marcus asked me to
## confirm the claim "FACTS is not designed to run to 2300" -- I could not, and it is false:
## EIGHT experiment configs in this very repo set `pyear_end: 2300` (the seven van Vuuren
## `src.*` arms and `fittedismip.ssp585`), with FittedISMIP / deconto21 / bamber19 / emulandice
## module sets, and SLEIP runs FACTS to 2300 as one of its seven emulators. The 2150 here is
## the horizon OUR base shared-climate arms use and the one the shipped comparison column
## reports at 2150; the *2300 twin keys below carry the SSP controls and van Vuuren scenarios
## to 2300 (climate-driven modules only past 2100).
## ⚠ A capability constant goes stale (memory `runnable_is_not_undrivable`): I read this
## constant's old comment -- "FACTS reaches 2150 at best" -- and reported it as a property of
## FACTS in four documents. Say RUN / RUNNABLE / IMPOSSIBLE, and check the machine.
## emulandice's 2100 cap IS a real module limit; it runs on the ten base keys at 2100 only.
FACTS_REPORT_MAX = 2150            # OUR shared-climate convention; NOT a FACTS ceiling
OUR_END_YEAR = 2300                # standing, memory `end_year_2300`
OHC_UNIT_TO_J = 1e22               # our cubes are ohc_1e22J


def _arg(flag, default):
    return next((a[len(flag):] for a in sys.argv[1:] if a.startswith(flag)), default)


def load_cube(kind, key):
    f = os.path.join(OBS, "fair_cube_%s_%s_raw.csv" % (kind, CUBE_KEY.get(key, key)))
    if not os.path.exists(f):
        raise SystemExit("no %s cube for %s at %s" % (kind, key, f))
    d = pd.read_csv(f)
    return d.year.values.astype(int), d.drop(columns=["year"])


def load_mean(kind, key):
    col = {"gmst": "gmst_C", "ohc": "ohc_1e22J"}[kind]
    d = pd.read_csv(os.path.join(OBS, "fair_mean_%s_%s.csv" % (kind, CUBE_KEY.get(key, key))))
    return dict(zip(d.year.astype(int), d[col].astype(float)))


def splice(years, raw, meanmap):
    """The Julia `convention()` for FORCING='spliced', per config column.

    raw: (nyear, ncfg) ndarray on `years`. Returns the same shape, spliced."""
    m = np.array([meanmap[y] for y in years], dtype=float)
    iref = np.flatnonzero((years >= REF[0]) & (years <= REF[1]))
    mref = m[iref].mean()
    cref = raw[iref, :].mean(axis=0)                       # per config
    out = np.where((years <= SPLICE_YEAR)[:, None],
                   m[:, None] * np.ones((1, raw.shape[1])),
                   mref + (raw - cref[None, :]))
    return out


# --- the config axis must be ONE axis across every scenario -----------------
def config_axis():
    """⚠ 'the same 841 configs' is a CLAIM ABOUT COLUMN ORDER, and it is checked, not
    assumed. If two cubes carry the configs in different orders then an index subsample
    silently selects a DIFFERENT climate per scenario, and the shared convention is false
    while every filename still says it is true."""
    ref_key, ref_cols = None, None
    for key, _f, _w, _l in SCENARIOS:
        for kind in ("gmst", "ohc"):
            _y, c = load_cube(kind, key)
            cols = list(c.columns)
            if ref_cols is None:
                ref_key, ref_cols = "%s/%s" % (key, kind), cols
            elif cols != ref_cols:
                d = [i for i, (a, b) in enumerate(zip(cols, ref_cols)) if a != b]
                raise SystemExit(
                    "[CONFIG-AXIS] %s/%s does not carry the same config order as %s "
                    "(%d column(s) differ, first at %d: %r vs %r). The shared-convention "
                    "claim would be FALSE. Refusing to build."
                    % (key, kind, ref_key, len(d), d[0], cols[d[0]], ref_cols[d[0]]))
    print("[CONFIG-AXIS] all %d scenarios x 2 fields carry an identical %d-config order"
          % (len(SCENARIOS), len(ref_cols)))
    return ref_cols


def build(key, nsamp, cfg_idx, outdir=None):
    facts_scen, lws_ssp, label = SCEN_BY_KEY[key]
    yr_g, cg = load_cube("gmst", key)
    yr_o, co = load_cube("ohc", key)
    if not np.array_equal(yr_g, yr_o):
        raise SystemExit("[YEARS] %s: gmst and ohc cubes on different year axes" % key)
    if yr_g.max() != OUR_END_YEAR:
        raise SystemExit("[YEARS] %s: cube ends at %d, expected the standing end year %d"
                         % (key, yr_g.max(), OUR_END_YEAR))
    ## [HORIZON] the 2300->2500 pad is a constant hold. It is only safe while nothing
    ## FACTS reports lives in it. Checked here so the assumption cannot rot silently.
    if horizon(key) > OUR_END_YEAR:
        raise SystemExit("[HORIZON] %s reports out to %d but our climate is real only to "
                         "%d -- the constant-hold pad would enter a reported year."
                         % (key, horizon(key), OUR_END_YEAR))

    G = splice(yr_g, cg.values.astype(float), load_mean("gmst", key))
    O = splice(yr_o, co.values.astype(float), load_mean("ohc", key))

    # [SPLICE-MATCH]: the same gate the Julia runs, against the same committed cube.
    ## CUBE_KEY, not key: the *2300 twins share their base key's committed cube, so the gate
    ## applies to ssp5852300 exactly as it does to ssp585 (it was silently "not applicable"
    ## on the twins before 2026-09-11).
    ref_f = os.path.join(OBS, "fair_cube_gmst_%s_spliced.csv" % CUBE_KEY.get(key, key))
    if os.path.exists(ref_f):
        R = pd.read_csv(ref_f)
        w = np.abs(G - R[list(cg.columns)].values.astype(float)).max()
        ok = w <= SPLICE_MATCH_TOL
        print("[SPLICE-MATCH] %-7s max |python splice - committed cube| = %.3e degC "
              "(tol %.0e)  %s" % (key, w, SPLICE_MATCH_TOL, "PASS" if ok else "FAIL"))
        if not ok:
            raise SystemExit("[SPLICE-MATCH] %s FAILED -- this is NOT the convention "
                             "Ladrillo and BRICK 2.0 are on, so it is not shared." % key)
    else:
        print("[SPLICE-MATCH] %-7s no committed python cube -- gate not applicable "
              "(it exists only for ssp585)" % key)

    # [BASELINE]: FACTS expects 1850-1900-centred SAT. Tolerance derived from the
    # ensemble's own spread over that window, not picked.
    ib = np.flatnonzero((yr_g >= 1850) & (yr_g <= 1900))
    base = G[ib, :].mean(axis=0)
    tol = 0.1 * G[ib, :].std(axis=0).mean()
    if np.abs(base).max() > max(tol, 1e-6):
        raise SystemExit("[BASELINE] %s: |1850-1900 mean| up to %.4f degC exceeds %.4f "
                         "(a tenth of the window's own sd). FACTS assumes 1850-1900 "
                         "centring." % (key, np.abs(base).max(), tol))
    print("[BASELINE]     %-7s |1850-1900 mean| max %.2e degC (tol %.2e)"
          % (key, np.abs(base).max(), max(tol, 1e-6)))

    G, O = G[:, cfg_idx], O[:, cfg_idx] * OHC_UNIT_TO_J
    ns = len(cfg_idx)

    def to_axis(a):
        out = np.empty((len(FACTS_YEARS), a.shape[1]), dtype=np.float32)
        for j in range(a.shape[1]):
            out[:, j] = np.interp(FACTS_YEARS, yr_g, a[:, j],
                                  left=a[0, j], right=a[-1, j])
        return out

    SAT, OHC = to_axis(G), to_axis(O)
    DOT = np.zeros_like(SAT)
    ny = len(FACTS_YEARS)
    outdir = outdir or os.path.join(FACTS_REPO, "experiments",
                                    "global.shared.%s.n%d" % (key, ns))
    os.makedirs(os.path.join(outdir, "input"), exist_ok=True)
    outdir = os.path.join(outdir, "input")
    prefix = "shared_%s" % key
    attrs = dict(
        Source="FaIR-2.2.4-calib1.6.0-CMIP7-injected",
        Scenario=facts_scen,
        Description=("%s; SHARED climate-driver convention: the same 841-config FaIR cube, "
                     "2014 splice, 1995-2014 reference that drive Ladrillo and BRICK 2.0. "
                     "cube key %s, FACTS scenario key %s, LWS socioeconomics %s, "
                     "n=%d of 841 configs."
                     % (label, key, facts_scen, lws_ssp, ns)),
        Centered="1850-1900 mean", Convention="spliced@2014, ref 1995-2014",
        CubeKey=key, LWSScenario=lws_ssp, FaIR_version="2.2.4 calib 1.6.0",
        PyearEnd=horizon(key))

    p = os.path.join(outdir, "%s_climate.nc" % prefix)
    with nc4.Dataset(p, "w", format="NETCDF4") as d:
        for k, v in attrs.items():
            setattr(d, k, v)
        d.createDimension("year", ny)
        d.createVariable("year", "i8", ("year",))[:] = FACTS_YEARS
        g = d.createGroup(facts_scen)
        g.createDimension("years", ny)
        g.createDimension("samples", ns)
        g.createVariable("years", "i8", ("years",))[:] = FACTS_YEARS
        g.createVariable("samples", "i8", ("samples",))[:] = np.arange(ns)
        for name, data, unit in [("surface_temperature", SAT, "degC"),
                                 ("ocean_heat_content", OHC, "J"),
                                 ("deep_ocean_temperature", DOT, "degC")]:
            vv = g.createVariable(name, "f4", ("years", "samples"), fill_value=np.nan)
            vv[:, :] = data
            vv.units = unit

    def write_flat(path, varname, data, unit):
        with nc4.Dataset(path, "w", format="NETCDF4") as d:
            for k, v in attrs.items():
                setattr(d, k, v)
            d.createDimension("samples", ns)
            d.createDimension("years", ny)
            d.createDimension("locations", 1)
            d.createVariable("samples", "i8", ("samples",))[:] = np.arange(ns)
            d.createVariable("years", "i8", ("years",))[:] = FACTS_YEARS
            d.createVariable("locations", "i8", ("locations",))[:] = [-1]
            d.createVariable("lat", "f8", ("locations",))[:] = [np.inf]
            d.createVariable("lon", "f8", ("locations",))[:] = [np.inf]
            vv = d.createVariable(varname, "f4", ("samples", "years", "locations"),
                                  fill_value=np.nan)
            vv[:, :, 0] = data.T
            vv.units = unit

    write_flat(os.path.join(outdir, "%s_gsat.nc" % prefix), "surface_temperature", SAT, "degC")
    write_flat(os.path.join(outdir, "%s_ohc.nc" % prefix), "ocean_heat_content", OHC, "J")

    def q(a, y, s=1.0):
        return np.percentile(a[np.flatnonzero(FACTS_YEARS == y)[0], :], [5, 50, 95]) / s
    s21, o21 = q(SAT, 2100), q(OHC, 2100, 1e21)
    print("  wrote %-46s n=%d  GMST@2100 %.2f [%.2f,%.2f] degC  OHC@2100 %.0f ZJ"
          % (os.path.relpath(outdir, FACTS_REPO), ns, s21[1], s21[0], s21[2], o21[1]))
    return dict(key=key, facts_scen=facts_scen, n=ns, gmst2100=s21[1], ohc2100=o21[1])


if __name__ == "__main__":
    NS = int(_arg("--nsamp=", "200"))
    cols = config_axis()
    ## ONE subsample, computed once, used for every scenario -- otherwise "the same
    ## configs" is false across the set even when each file is individually correct.
    CFG_IDX = np.unique(np.linspace(0, len(cols) - 1, NS).round().astype(int))
    print("[SUBSAMPLE]    n=%d of %d configs, evenly spaced, IDENTICAL across all scenarios "
          "(first 3: %s)" % (len(CFG_IDX), len(cols), [cols[i] for i in CFG_IDX[:3]]))
    keys = ([k for k, _f, _w, _l in SCENARIOS] if "--all" in sys.argv[1:]
            else [_arg("--scen=", "ssp245")])
    rows = [build(k, NS, CFG_IDX) for k in keys]
    print("\n%-7s %-9s %6s %12s %10s" % ("key", "scen key", "n", "GMST@2100", "OHC@2100"))
    for r in rows:
        print("%-7s %-9s %6d %12.2f %10.0f"
              % (r["key"], r["facts_scen"], r["n"], r["gmst2100"], r["ohc2100"]))
    print("\nDriving path: scenario emissions -> FaIR 2.2.4 calib 1.6.0 (841 cfgs) -> GMST+OHC "
          "-> FACTS. No module on it reads climate from a scenario label; only landwaterstorage "
          "(population pathway, LWSScenario) and emulandice (2100, per-SSP-trained) consult one.")
