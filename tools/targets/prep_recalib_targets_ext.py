#!/usr/bin/env python3
"""
prep_recalib_targets_ext.py -- build Ladrillo's calibration targets,
data/observations/calibration_targets.csv (research name
outputs/recalib_targets_ext.csv).

PORT of the research `python/prep_recalib_targets_ext.py` with no flags
(SLR-RFF-BRICK @ bc50130; the L27 target, Frederikse AIS). Every line that
computes a target column is verbatim. What changed:
  * paths come from --sources (data/sources/) and --work (the build directory)
    instead of the research tree; the research raw/ and data/observations/
    file names are kept under data/sources/<provider>/;
  * the IMBIE-2026 AIS branch (--ais-imbie2026, the closed L28-L35 arms) is
    dropped, so the script builds the L27 target only;
  * a missing GRACE land-water extension is now an ERROR. The research script
    printed a warning and fell back to holding LWS flat from 2018, which builds
    a different target without failing;
  * dropped because nothing downstream reads them: the IMBIE 2023 cross-check
    and its two input files, the splice diagnostic PNG, the
    recalib_targets_ext_sources.csv sidecar;
  * the provenance sidecar records the Ladrillo commit and the sha256 of every
    input instead of the research commit.

Inputs (under --sources): frederikse2020/ (workbook, 5000-member ensemble),
dangendorf2024/ (annual GMSL extraction, the v2 global file), grace/ (the JPL
Antarctic and Greenland series), glambie/calendar_years/ (global and region
5), noaa/ (thermosteric 0-2000 m, STAR altimetry annual means); under --work:
outputs/lws_grace_extension_L14.csv (build_lws_grace_extension.py).

The research docstring follows, unchanged.

Prep EXTENDED recalibration targets: Frederikse 2020 (1900-2018) spliced with
modern reconciled products to extend each SLR component past 2018.

Motivation (Marcus, 2026-06-13): test how post-2018 observations -- especially
the post-2020 Antarctic GRACE-FO "pause" -- shift the BRICK-Mengel calibration,
and extend Greenland / thermal-expansion / glaciers at the same time.

Extension data (all reconciled multi-method products; see README in raw/):
  - AIS : ⭐ SINCE 2026-09-21 (Marcus): IMBIE 2026 (Otosaka et al., Sci Data 13:1301; PDC 10.5285/128c5e33)
          reconciled Antarctic mass balance 1979-2023 REPLACES Frederikse's AIS over those years; Frederikse
          1900-1978 is offset-matched onto it over IMBIE_JOIN_WIN; GRACE-FO JPL mascon RL06.3Mv4 (DOI
          10.5067/TEMSC-3JC634) is offset-matched onto IMBIE over 2003-2023 and fills 2024-2026 only.
          `--ais-frederikse` restores the pre-09-21 target (Frederikse to 2018 + GRACE from 2019), byte-identical.
          Why: diag_imbie2026_vs_targets.py -- the old target sat 22 % below the reconciled record on the 1992-2020
          rate (-1.3 sigma with both bars), every window since 1979 low; CHANGELOG 09-21f/g.
  - GIS : GRACE-FO JPL mascon RL06.3Mv4 Greenland mass (2002-2026), same DOI
  - GSIC: GlaMBIE 2025 glacier mass (2000-2023), DOI 10.5904/wgms-glambie-2024-07,
          SCOPE-MATCHED = global MINUS region 5 (r5 is in the GIS target), region 19 KEPT
          (2026-08-06 decision; see the GSIC block below + memo_2026-08-05 §2c/§3)
  - TE  : NOAA NCEI World-Ocean 0-2000m thermosteric sea level (2005-2025)
  - TOT : real Dangendorf 2024 GMSL reconstruction (1900-2021) spliced with NOAA STAR
          altimetry (2022-2024). [M3 rework 2026-07-20: was Frederikse-total + STAR;
          Frederikse is now COMPONENTS-only, the total is the independent Dangendorf.]
  - IMBIE 2023 (Otosaka et al., 1992-2020) AIS+GIS used only as an INDEPENDENT
    cross-check on the GRACE splices over the overlap (not fed to the fit).

Component band σ (2026-07-20): from the Frederikse 5000-member weighted ensemble
(FRED_ENS_NC), re-referenced per member -- the statistically correct band that shrinks
toward the 1995-2005 window, finishing the 2026-07-19 σ fix.

Method (each component, matching prep_recalib_targets.py conventions):
  1. Frederikse component re-referenced to the common 1995-2005 window mean (cm).
  2. Modern product converted to cm sea-level-equivalent (mass: -Gt/GT_PER_CM_SLE;
     steric/altimetry: mm/10), then OFFSET-matched (pure level shift, no rescale --
     both measure the same physical SLE) so its mean over the component's overlap
     window equals Frederikse's mean over the same window.
  3. Frederikse years 1900-2018, then the spliced modern product for 2019..end.

Decisions baked in (Marcus 2026-06-13, see handoff_2026-06-13_..._ais_extension):
  - Total extended with NOAA STAR altimetry (budget stays closed).
  - LWS 2019-2023 is REAL GRACE/GRACE-FO data (build_lws_grace_extension.py), wired in
    2026-08-25 alongside the A6 amp re-centring; held flat only past the last real year.
    ⚠ its BAND WIDTH is still the 2018 half-width -- a real centre with a held bar. FLAGGED.
  - the budget-closure sigma is TREND-EXTENDED past the ensemble (2026-08-25), not held
    flat: the hold was measurably anti-conservative (+0.0216 cm/yr, z=+6.08).
  - Heterogeneous end years per component (AIS/GIS 2026, GSIC 2023, TE 2025, TOT 2024);
    the per-series AR(1) likelihood in calibrate_mcmc_ext.jl handles different lengths.

Output: outputs/recalib_targets_ext.csv (1900-2026; NaN where a component has no data)
        outputs/recalib_targets_ext_splice_diagnostic.png
"""
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

import argparse
_ap = argparse.ArgumentParser(description="Ladrillo calibration targets (L27)")
_ap.add_argument("--sources", required=True, help="data/sources directory")
_ap.add_argument("--work", required=True, help="build directory (reads outputs/lws_grace_extension_L14.csv, "
                                               "writes outputs/recalib_targets_ext.csv)")
_ap.add_argument("--commit", default="unknown", help="Ladrillo commit, for the provenance sidecar")
_ARGS = _ap.parse_args()
AIS_SOURCE = "frederikse"     # PORT: the only build Ladrillo ships (the research --ais-imbie2026 arm is dropped)

# ⚠⚠ WHICH POSTERIOR DOES THIS TARGET REPRODUCE? The default was IMBIE 2026 from 2026-09-21
# (commit c247e06) to 2026-09-30, while the shipped posterior stayed on Frederikse -- so a bare
# rebuild did NOT reproduce the champion. Marcus 2026-09-30: L27 is the default everywhere, so the
# default is back on Frederikse; the IMBIE build needs --ais-imbie2026 (its last output is kept as
# outputs/recalib_targets_ext_imbie2026*). Say which, out loud, before anything is written.
AIS_SOURCE_REPRODUCES = {
    "frederikse": "L27 and earlier (the SHIPPED champion and the GMD draft's posterior)",
    "imbie2026":  "L28-L33 (the IMBIE-2026 arms) -- NOT the shipped champion",
}
print("=" * 100)
print(f"  AIS CALIBRATION TARGET = {AIS_SOURCE.upper()}"
      f"{'  (NOT the shipped L27 target -- that is the default, no flag)' if AIS_SOURCE != 'frederikse' else ''}")
print(f"  reproduces: {AIS_SOURCE_REPRODUCES[AIS_SOURCE]}")
print("=" * 100)

SRC = _ARGS.sources
OUT_CSV = os.path.join(_ARGS.work, "outputs/recalib_targets_ext.csv")
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

# Frederikse 5000-member weighted component ENSEMBLE (redistributed in Dangendorf's Zenodo
# record 10621070). This is the object the 2026-07-19 σ-fix flagged as missing: it lets us
# compute the statistically correct re-referenced per-year band σ = weighted sd of members
# each re-referenced to the BASE window, which SHRINKS toward the window (reconstruction
# errors are correlated in time) instead of the raw-band-width over-statement.
FRED_ENS_NC = os.path.join(SRC, "frederikse2020/frederikse2020_GMSL_ensembles.nc")
# ensemble var -> target key
FRED_ENS_VAR = {"AIS": "ais", "Glaciers": "gsic", "GrIS": "gis",
                "Steric": "steric", "TWS": "lws", "GMSL": "dang"}

# ---- budget-closure inflation of the TOTAL target's sigma (Marcus, 2026-08-12) ----
# Gate 3.1 established that the five component targets sum +0.74 cm above the independent
# total target over 1950-1980, and that this is Frederikse 2020's OWN mid-century budget
# non-closure (our closure sits at z = +0.006 against their 5000-member ensemble median).
# RULING: carry it as uncertainty on the total, across the WHOLE span rather than a chosen
# window, using the ensemble's OWN per-year closure spread -- no decay function is fitted
# and no window edge is chosen. Members are re-referenced to BASE exactly as the targets
# are, so sigma is smallest at the anchor and grows away from it in both directions; that
# is the correct structure for a level anomaly scored in that frame, and it is why the
# inflation is 1.37x at 1900, 1.11x at 2000 and 1.90x at 2018 rather than monotone.
# CAVEAT ON RECORD: the total channel also carries a sampled AR(1) with rho ~ 0.97, so an
# anchor-shaped sigma may partly double-count level correlation. Flagged, not corrected.
CLOSURE_COMPONENT_VARS = ["AIS", "Glaciers", "GrIS", "Steric", "TWS"]
CLOSURE_TOTAL_VAR = "GMSL"
CLOSURE_SIG_COL = "dang_closure_sig"
# The ensemble ends in 2018 and the total target runs to 2024. Hold flat past the ensemble,
# the same FLAGGED convention already used for LWS.
# ⚠ FLIPPED FALSE 2026-08-25: the hold is measurably anti-conservative (see the block where
# it is used). True reproduces every pre-L15 target file exactly.
CLOSURE_HOLD_FLAT_PAST_ENSEMBLE = False
# The window the closure-sigma trend is fitted over -- the last decade of the ensemble.
CLOSURE_TREND_WIN = (2009, 2018)
# The tag of the GRACE LWS extension to wire in. Its builder is
# python/build_lws_grace_extension.py; it does NOT depend on a posterior, the tag is only
# the filename convention it was written under.
LWS_GRACE_TAG = "L14"

# ---- constants ----
BASE_Y0, BASE_Y1 = 1995, 2005          # common re-reference window (11 yr)
FIT_Y0           = 1900                 # fit start
EXT_Y1           = 2026                 # widest extension (AIS/GIS GRACE end)
GT_PER_CM_SLE    = 3620.0              # 362 Gt/mm * 10 (ocean-area based; <1% on increments)
ALT_SIGMA_MM     = 4.0                 # nominal altimetry annual GMSL 1-sigma (NOAA STAR sigma col empty)
# ⭐ RULED 2026-09-09c (Marcus): ALT_SIGMA_MM STAYS AT 4.0. Measured against the across-product
# spread IN THE SPLICE YEARS (2022-24) it is 1.9-3.2 mm, and Dangendorf's own SE at the join is
# 2.68 -- so 4.0 is CONSERVATIVE and replacing it would TIGHTEN the target, the opposite of what
# the 09-09 scoping note implied. The scope's "1.4-11x" was the ERA-WIDE range, evaluated over
# the wrong window ([[like_for_like_forcing]] applies to a sigma's window too).
# Evidence: python/diag_modern_splice_altimetry.py; notes/handoff_2026-09-09b... section 4.
# per-component overlap window for the offset-match splice
OVERLAP = {"ais": (2003, 2018), "gis": (2003, 2018), "gsic": (2003, 2018),
           "steric": (2005, 2018), "dang": (2003, 2018)}
FRED = os.path.join(SRC, "frederikse2020/frederikse2020_global_basin_timeseries.xlsx")
FRED_MAP = {"Antarctic Ice Sheet": "ais", "Glaciers": "gsic",
            "Greenland Ice Sheet": "gis", "Steric": "steric",
            "Terrestrial Water Storage": "lws"}


def reref(s, win):
    return s - s.loc[win[0]:win[1]].mean()


def _wquantile_free(vals, w, q):
    idx = np.argsort(vals)
    xs, ws = vals[idx], w[idx]
    cw = (np.cumsum(ws) - 0.5 * ws) / ws.sum()
    return np.interp(q, cw, xs)


def load_ensemble_sigma():
    """Per-year band σ (cm) for each component from the Frederikse 5000-member ensemble.

    σ_t = weighted sd across members of (x_{m,t} − mean_{window} x_{m,·}), i.e. each member
    re-referenced to BASE_Y0..BASE_Y1 with its OWN single offset (band width preserved, but
    the spread correctly collapses toward the reference window). Returns {tgt: Series(cm)}
    indexed by ensemble years (1900-2018).
    """
    ds = xr.open_dataset(FRED_ENS_NC)
    ey = ds["time"].values.astype(int)
    w = ds["likelihood"].values.astype(float)
    win = (ey >= BASE_Y0) & (ey <= BASE_Y1)
    sig = {}
    for var, tgt in FRED_ENS_VAR.items():
        m = ds[var].values                                   # (member, year), mm
        m = m - m[:, win].mean(axis=1, keepdims=True)        # per-member single offset
        mean = np.average(m, axis=0, weights=w)
        sd = np.sqrt(np.average((m - mean) ** 2, axis=0, weights=w))
        sig[tgt] = pd.Series(sd / 10.0, index=ey)            # mm -> cm
    return sig


def load_closure_sigma():
    """Per-year weighted sd (cm) of Frederikse's OWN budget closure,
    (sum of the five components) - (their GMSL), across the 5000 members.

    Each member is re-referenced to BASE_Y0..BASE_Y1 with its own single offset,
    identically to load_ensemble_sigma() and to the targets themselves, so this
    sigma is in the same frame as the residual the likelihood scores. Indexed by
    ensemble year (1900-2018)."""
    ds = xr.open_dataset(FRED_ENS_NC)
    ey = ds["time"].values.astype(int)
    w = ds["likelihood"].values.astype(float)
    win = (ey >= BASE_Y0) & (ey <= BASE_Y1)

    def reref_members(var):
        m = ds[var].values / 10.0                            # mm -> cm
        return m - m[:, win].mean(axis=1, keepdims=True)

    close = (sum(reref_members(v) for v in CLOSURE_COMPONENT_VARS)
             - reref_members(CLOSURE_TOTAL_VAR))             # (member, year)
    mean = np.average(close, axis=0, weights=w)
    sd = np.sqrt(np.average((close - mean) ** 2, axis=0, weights=w))
    return pd.Series(sd, index=ey)


def annual_mean(t, v, ylo, yhi):
    """Calendar-year mean of an (irregular) monthly series t (decimal yr) -> v."""
    yr = np.floor(t).astype(int)
    out = {}
    for y in range(ylo, yhi + 1):
        m = yr == y
        if m.sum() >= 6:                # require >=6 months for a usable annual mean
            out[y] = float(np.mean(v[m]))
    return pd.Series(out)


def splice_offset(modern_cm, fred_cm, win):
    """Level shift so modern mean == Frederikse mean over the overlap window."""
    yrs = [y for y in range(win[0], win[1] + 1) if y in modern_cm.index and y in fred_cm.index]
    return float(fred_cm.loc[yrs].mean() - modern_cm.loc[yrs].mean()), yrs


# ============================================================ Frederikse (rel window)
g = pd.read_excel(FRED, "Global").rename(columns={"Unnamed: 0": "year"}).set_index("year")
years = np.arange(FIT_Y0, EXT_Y1 + 1)
out = pd.DataFrame({"year": years}).set_index("year")
ens_sig = load_ensemble_sigma()            # per-year band σ (cm), 1900-2018, per component
fred = {}
for fname, tgt in FRED_MAP.items():
    # VALUE = Frederikse component mean (re-referenced with the mean's single offset).
    off = (g[f"{fname} [mean]"] / 10.0).loc[BASE_Y0:BASE_Y1].mean()    # mm->cm, one offset
    mean_s = (g[f"{fname} [mean]"] / 10.0 - off).reindex(years)
    out[tgt] = mean_s.values
    fred[tgt] = (g[f"{fname} [mean]"] / 10.0) - off
    # BAND σ (2026-07-20): use the ENSEMBLE per-year sd, not the mean/lower/upper band width.
    # This finishes the 2026-07-19 σ fix. The old raw-band approach: (i) needed the mean's
    # single offset to avoid inverting lo>hi (the shrink-by-constant bug), and even fixed it
    # OVER-stated σ near the reference window because the raw band doesn't collapse there.
    # sd over the re-referenced ensemble members DOES collapse toward 1995-2005 (correlated
    # reconstruction errors) -- the statistically correct quantity. _lo/_hi are written as
    # value ∓ 1.645·σ so calibrate_mcmc_ext.jl's ϵband=(hi-lo)/(2·1.645) recovers σ exactly.
    sig = ens_sig[tgt].reindex(years)
    out[tgt + "_lo"] = (mean_s - 1.645 * sig).values
    out[tgt + "_hi"] = (mean_s + 1.645 * sig).values

# "Total" term (rel window). M3 REWORK 2026-07-20 (Marcus): the total is now the REAL
# Dangendorf 2024 reconstruction (1900-2021), NOT Frederikse (which the old file secretly
# was). Dangendorf is an independent-of-Frederikse 20th-c. tide-gauge reconstruction; the
# tension diagnostic (diag_dangendorf_vs_frederikse.py) shows it sits inside Frederikse's
# 5-95% at every trend window (mid-century 6.8th pctl) and agrees with altimetry better in
# the satellite era -- so keeping Frederikse COMPONENTS + a Dangendorf TOTAL introduces <1σ
# of trend inconsistency. σ (2026-08-07, Marcus): Dangendorf's OWN per-year SE from the
# CORRECTED Global_v2.nc (S. Dangendorf pers.comm. 2026-08-07; the original Zenodo Global.nc
# was slot-shifted -- GMSL slots held barystatic -- which is why the 2026-07-20 rework had
# substituted the Frederikse ensemble GMSL sd; that workaround claimed conservatism, but the
# corrected SE is 1.3-2x LARGER over 1900-2010 and smaller post-2015, so the substitution is
# RETIRED). SE units are METERS in the file ("the error ... is 3 mm in 2021 ... given in m").
# Covers 1900-2021, so no hold-flat extrapolation; STAR years (2022+) get ALT_SIGMA_MM below.
DANG_V2_NC = os.path.join(SRC, "dangendorf2024/dangendorf2024_KalmanSmootherHR_Global_v2.nc")
d = pd.read_csv(os.path.join(SRC, "dangendorf2024/dangendorf2024_gmsl_annual.csv")).set_index("year")
dang_val = reref(d["gmsl_mm"] / 10.0, (BASE_Y0, BASE_Y1))                 # cm, 1900-2021
out["dang"]  = dang_val.reindex(years).values
fred["dang"] = dang_val
with xr.open_dataset(DANG_V2_NC) as dv2:
    dang_sig = pd.Series(dv2["GMSLHRSE"].values.ravel() * 100.0,          # m -> cm
                         index=dv2["t"].values.ravel().astype(int)).reindex(years)
assert abs(dang_sig.loc[2021] - 0.268) < 0.01, "v2 SE(2021) should be ~2.68 mm (units check)"
out["dang_sig"] = dang_sig.values
out["dang_lo"]  = (dang_val.reindex(years) - 1.645 * dang_sig).values
out["dang_hi"]  = (dang_val.reindex(years) + 1.645 * dang_sig).values

# budget-closure inflation of the total's sigma (see the constant block; Marcus 2026-08-12).
# Emitted as its OWN column: dang_lo/dang_hi stay the Dangendorf band, and
# calibrate_mcmc_ext.jl adds this in quadrature alongside the LWS band term.
closure_sig = load_closure_sigma()
_ens_last = int(closure_sig.index.max())
closure_col = closure_sig.reindex(years)
if CLOSURE_HOLD_FLAT_PAST_ENSEMBLE:
    closure_col.loc[_ens_last + 1:] = closure_sig.loc[_ens_last]
else:
    # ⚠ MEASURED 2026-08-25: HOLDING THIS FLAT IS ANTI-CONSERVATIVE. Over the last decade of
    # the ensemble the closure spread GROWS at +0.0216 cm/yr, and that slope is RESOLVED --
    # z = +6.08 against an AR(1)-inflated bar. Freezing it at the 2018 value therefore makes
    # the TOTAL target's error bar too TIGHT in exactly the years where we have least
    # information about budget closure, over-weighting the recent total in the fit.
    #
    # Extending an UNCERTAINTY on a resolved trend is not the same move as extrapolating a
    # LEVEL, which `postsplice_halving_priced` records failing here (two parametric arms
    # agreed with each other and both missed the sign). Erring WIDE past the end of the data
    # is the safe direction for an error bar; erring wide on a level is not.
    #
    # ⚠ AND THE HOLD IS NOT FAILING IN PRACTICE, which is why this is a WARN-level fix and
    # not a correction: the realised budget residual over 2019-2023, computed from the modern
    # extension products with the GRACE LWS, has max |residual| 0.532 cm = 0.69x the held
    # 0.775. Its MEAN shifts (-0.375 cm against +0.119 pre-2019) while its sd SHRINKS
    # (0.123 against 0.218) -- a systematic splice offset, not extra scatter. So the widened
    # bar is insurance, not a repair.
    _fit = closure_sig.loc[CLOSURE_TREND_WIN[0]:CLOSURE_TREND_WIN[1]]
    _b, _a = np.polyfit(_fit.index.values.astype(float), _fit.values, 1)
    _fut = np.arange(_ens_last + 1, years[-1] + 1, dtype=float)
    closure_col.loc[_ens_last + 1:] = np.maximum(_a + _b * _fut, closure_sig.loc[_ens_last])
    print(f"\nclosure sigma EXTENDED on the {CLOSURE_TREND_WIN[0]}-{CLOSURE_TREND_WIN[1]} "
          f"trend {_b:+.5f} cm/yr (z=+6.08): {closure_sig.loc[_ens_last]:.4f} cm at "
          f"{_ens_last} -> {closure_col.loc[years[-1]]:.4f} cm at {years[-1]} "
          f"(floored at the last ensemble value)")
out[CLOSURE_SIG_COL] = closure_col.values

# ============================================================ modern products -> cm SLE
modern = {}      # tgt -> (Series cm, Series sigma_cm)

# --- AIS / GIS : GRACE-FO mascon (Gt anomaly) ---
for tgt, fn in [("ais", "grace/grace_antarctica_mass.txt"), ("gis", "grace/grace_greenland_mass.txt")]:
    raw = pd.read_csv(os.path.join(SRC, fn), sep=r"\s+", comment="H", header=None,
                      names=["t", "m", "s"]).dropna()
    sle = annual_mean(raw.t.values, -raw.m.values / GT_PER_CM_SLE, 2002, EXT_Y1)
    sig = annual_mean(raw.t.values,  raw.s.values / GT_PER_CM_SLE, 2002, EXT_Y1)
    modern[tgt] = (sle, sig)

# --- GSIC : GlaMBIE glacier mass, SCOPE-MATCHED (annual Gt change -> cumulative) ---
# 2026-08-06 scope decision (Marcus; memo_2026-08-05 §3 choice 1c): the GSIC component owns
# RGI regions 1-18-minus-5 PLUS region 19. Region 5 (Greenland periphery) lives in the GIS
# target (Frederikse GrIS = Kjeldsen/Mouginot + r5; GRACE mascon Greenland includes it), so
# the GlaMBIE GLOBAL series must have r5 SUBTRACTED or the splice double-counts it against
# GIS. Region 19 is RETAINED (inventory convention V = 0.290 m SLE). NB the Frederikse
# 1900-2018 glacier segment lacks r19 flow (documented zero in his Methods) -- small known
# bias (~0.05 mm/yr modern, less earlier). ⚠ r19 is a zero in Frederikse's GLACIER column
# only, NOT "everywhere else": every gravimetric AIS source (his GRACE era, the GRACE-FO
# splice, IMBIE 2026's gravimetry) includes the connected Antarctic periphery, by a fraction
# no product reports. Ceiling on the resulting double count = all of GlaMBIE r19 over the
# gravimetry years, 0.075 cm 2002-2021 -- deliverables/obs_consistency_vs_imbie2026.md §5.
# Error: r5's sigma quadrature-ADDED to the global sigma (conservative; r5 ~35 vs global
# ~100 Gt/yr, so the change is minor).
gl = pd.read_csv(os.path.join(SRC, "glambie/calendar_years/0_global.csv"))          # = research raw/glambie_global_glacier_mass.csv
r5 = pd.read_csv(os.path.join(SRC, "glambie/calendar_years/5_greenland_periphery.csv"))  # = raw/glambie_r5_greenland_periphery.csv
assert np.array_equal(gl.start_dates.values, r5.start_dates.values), "GlaMBIE global/r5 year mismatch"
gl["yr"] = np.floor(gl.start_dates).astype(int)
gl_net_gt = gl.combined_gt.values - r5.combined_gt.values         # global minus r5 (r19 kept)
gl_net_err = np.sqrt(gl.combined_gt_errors.values ** 2 + r5.combined_gt_errors.values ** 2)
gl_cum = pd.Series(-gl_net_gt, index=gl.yr.values).cumsum() / GT_PER_CM_SLE   # cm SLE, cumulative
gl_sig = pd.Series(np.sqrt((gl_net_err ** 2).cumsum()) / GT_PER_CM_SLE, index=gl.yr.values)
modern["gsic"] = (gl_cum, gl_sig)
print(f"GSIC splice scope: GlaMBIE global - r5 (r19 kept); 2000-2023 cumulative "
      f"{gl_cum.iloc[-1]:+.2f} cm (global would be {-gl.combined_gt.sum()/GT_PER_CM_SLE:+.2f})")

# --- TE : NOAA NCEI World-Ocean 0-2000m thermosteric SL (mm) ---
st = pd.read_csv(os.path.join(SRC, "noaa/noaa_thermosteric_w0-2000m_yearly.dat"),
                 sep=r"\s+", skiprows=1, header=None,
                 names=["t", "WO", "WOse", "NH", "NHse", "SH", "SHse"])
st["yr"] = np.floor(st.t).astype(int)
modern["steric"] = (pd.Series((st.WO / 10.0).values, index=st.yr.values),
                    pd.Series((st.WOse / 10.0).values, index=st.yr.values))

# --- TOTAL : NOAA STAR altimetry GMSL (mm) ---
# ⭐⭐ RULED 2026-09-09c (Marcus): AT THE NEXT TARGET REBUILD THIS ANCHOR BECOMES THE IGCC
# THREE-PRODUCT ALTIMETRY ENSEMBLE, not NOAA STAR. NOT APPLIED HERE, DELIBERATELY: the same
# ruling scoped the 09-09c update as PROVENANCE-ONLY, so no target was rebuilt and no chain
# was run. Editing this line without regenerating outputs/recalib_targets_ext.csv would leave
# the script disagreeing with its own output -- do not half-apply it.
#   why IGCC: STAR is a single product, and the one IGCC dropped (for AVAILABILITY, not
#   quality); offset-matched over 2003-2018 it sits 3.67 mm = 0.92 sigma BELOW the IGCC
#   ensemble at 2024, with every IGCC product above it in every spliced year. IGCC also
#   reaches 2025, one year past STAR.
#   reach: 3 years directly (2022-24), but the 2003-2018 offset match means the choice moves
#   those years twice over.
#   ⚠ the 3-product ensemble USES the `NOAA`-headed col 3 of altimetry_indiv_estimates.csv,
#   which a standing note says not to use until its header/product mismatch is settled. The
#   ruling lifts that hold for this purpose; measured cost of honouring it instead is 0.02 mm
#   at 2024 (2-product mean 85.90 vs 3-product 85.92), so nothing here turns on it.
# Evidence: python/diag_modern_splice_altimetry.py; notes/handoff_2026-09-09b... section 4-5.
alt = pd.read_csv(os.path.join(SRC, "noaa/nasa_gmsl_annual.csv")).set_index("year")
modern["dang"] = (alt["value"] / 10.0, pd.Series(ALT_SIGMA_MM / 10.0, index=alt.index))

# ============================================================ splice modern products into out
# Splice starts the year AFTER each historical series ends: components extend Frederikse
# (ends 2018) from 2019; the TOTAL extends the real Dangendorf reconstruction (ends 2021)
# with NOAA STAR altimetry from 2022 -- so STAR fills only the years Dangendorf lacks.
SPLICE_FROM = {"ais": 2019, "gis": 2019, "gsic": 2019, "steric": 2019, "dang": 2022}
print(f"{'comp':6s} {'overlap':12s} {'offset_cm':>9s} {'from':>5s} {'end_yr':>6s}  increment to end (cm)")
splices = {}
for tgt in ["ais", "gis", "gsic", "steric", "dang"]:
    mod_cm, mod_sig = modern[tgt]
    off, yrs = splice_offset(mod_cm, fred[tgt], OVERLAP[tgt])
    y0 = SPLICE_FROM[tgt]
    spl = mod_cm + off
    end = int(spl.index.max())
    splices[tgt] = (spl, mod_sig)
    for y in range(y0, end + 1):
        if y in spl.index:
            out.loc[y, tgt] = spl[y]
            s = mod_sig.get(y, np.nan)
            s = 0.05 if not np.isfinite(s) else max(s, 0.01)
            out.loc[y, tgt + "_lo"] = spl[y] - 1.645 * s
            out.loc[y, tgt + "_hi"] = spl[y] + 1.645 * s
    if tgt == "dang":
        for y in range(y0, end + 1):
            if y in spl.index:
                out.loc[y, "dang_sig"] = ALT_SIGMA_MM / 10.0
    inc = spl[end] - fred[tgt].get(y0 - 1, np.nan)
    print(f"{tgt:6s} {str(OVERLAP[tgt]):12s} {off:9.3f} {y0:5d} {end:6d}  {inc:+.3f}")

# LWS 2019+: REAL DATA (GRACE/GRACE-FO), replacing the hold-flat fiat -- WIRED IN 2026-08-25.
#
# `python/build_lws_grace_extension.py` builds LWS from the JPL mascons minus GlaMBIE
# glaciers. It was written 2026-08-24 and deliberately NOT wired in, because swapping a
# fitted target is a recalibration trigger and the level alone did not justify one (GRACE
# differs from the hold by mean +0.018 cm, sd 0.076, max 0.123 -- the hold got the LEVEL
# right and removed only the interannual variance). Marcus: "wait to recalibrate until we
# have something else worth recalibrating." The A6 amp re-centring is that something, so
# this rides along with it.
#
# ⚠ WHY IT MATTERS EVEN THOUGH THE LEVEL BARELY MOVES. LWS enters the likelihood TWICE on
# the total stream -- `calibrate_mcmc_ext.jl:1387` adds observed LWS to the MODELLED total,
# and `:568` folds its band into the total's observation sigma -- so for 2019+ the fit was
# being told, AS DATA WITH A REAL-DATA ERROR BAR, that land water contributed exactly its
# 2018 value every year.
#
# ⚠ THE CENTRE IS NOW REAL; THE BAND WIDTH IS STILL HELD, and that is a deliberate,
# smaller fiat. GRACE mascons carry their own uncertainty but the extension builder does not
# emit one, and the GRACE/Frederikse overlap is 2 years -- too few to fit a width from. The
# 2018 half-width is carried forward (TWS uncertainty does not shrink), so the fit now sees
# a REAL value with a HELD error bar instead of a held value with a real error bar. Flagged.
LWS_GRACE_CSV = os.path.join(_ARGS.work, "outputs", f"lws_grace_extension_{LWS_GRACE_TAG}.csv")
_lws_half = (out.loc[2018, "lws_hi"] - out.loc[2018, "lws_lo"]) / 2.0
_n_real = 0
# PORT: the research script fell back to holding LWS flat from 2018 when this file was missing
# (a warning, then a DIFFERENT target). A missing extension is now an error.
if not os.path.exists(LWS_GRACE_CSV):
    sys.exit(f"ABORT: {LWS_GRACE_CSV} not found -- run build_lws_grace_extension.py first "
             f"(the research script held LWS flat from 2018 here, silently building a different target)")
if True:
    _g = pd.read_csv(LWS_GRACE_CSV).set_index("year")
    # Only years the builder marks as REAL (is_held_flat False) are used; its own tail hold
    # is not imported, or we would be re-importing a fiat under a different name.
    _real = _g.loc[~_g.is_held_flat.astype(bool), "lws_grace_cm"].dropna()
    for y, v in _real.items():
        if 2019 <= y <= EXT_Y1:
            out.loc[y, "lws"] = v
            out.loc[y, "lws_lo"] = v - _lws_half
            out.loc[y, "lws_hi"] = v + _lws_half
            _n_real += 1
    _last_real = int(_real.index.max()) if len(_real) else 2018
# Past the last REAL GRACE year the hold-flat convention still applies, now anchored on the
# last real value rather than on 2018.
_vlast = out.loc[_last_real, "lws"]
for suf in ["", "_lo", "_hi"]:
    ref = out.loc[_last_real, "lws" + suf]
    out.loc[_last_real + 1:EXT_Y1, "lws" + suf] = ref
print(f"\nLWS: {_n_real} year(s) of REAL GRACE data wired in (2019-{_last_real}); "
      f"held flat {_last_real + 1}-{EXT_Y1} at {_vlast:.4f} cm; band half-width "
      f"{_lws_half:.4f} cm carried from 2018 (FLAGGED)")

# budget-closure inflation, reported against the FINAL sigma columns (post-splice, post
# LWS hold-flat) so the printed numbers are the ones calibrate_mcmc_ext.jl will build.
_eps_lws = ((out["lws_hi"] - out["lws_lo"]) / (2 * 1.645)).clip(lower=0.05)
_cur = np.sqrt(out["dang_sig"] ** 2 + _eps_lws ** 2)
_new = np.sqrt(_cur ** 2 + out[CLOSURE_SIG_COL] ** 2)
print(f"\nTotal-target sigma inflated by Frederikse's own budget-closure spread "
      f"(ensemble {closure_sig.index.min()}-{_ens_last}, "
      f"{'held flat' if CLOSURE_HOLD_FLAT_PAST_ENSEMBLE else 'trend-extended'} after):")
for y in [1900, 1950, 1980, 2000, 2018, 2024]:
    if y in out.index and np.isfinite(_cur.loc[y]) and np.isfinite(out.loc[y, "dang"]):
        print(f"  {y}  closure sd {out.loc[y, CLOSURE_SIG_COL]:5.3f} cm   "
              f"sigma {_cur.loc[y]:5.3f} -> {_new.loc[y]:5.3f} cm  "
              f"({_new.loc[y] / _cur.loc[y]:.2f}x)")

out.reset_index().to_csv(OUT_CSV, index=False)
print(f"\nWrote {OUT_CSV}  ({len(out)} yrs {FIT_Y0}-{EXT_Y1})")

# ============================================================ provenance sidecar (PORT)
# The research script also wrote recalib_targets_ext_sources.csv, an IMBIE 2023 cross-check and a
# splice PNG; nothing downstream reads them. The provenance stamp is kept, with the Ladrillo commit
# and the sha256 of every input (the research md5 line is replaced by sha256).
import datetime as dt, hashlib
def _sha(p):
    with open(p, "rb") as _fh:
        return hashlib.sha256(_fh.read()).hexdigest()
_inputs = [FRED, FRED_ENS_NC, DANG_V2_NC, os.path.join(SRC, "dangendorf2024/dangendorf2024_gmsl_annual.csv"),
           os.path.join(SRC, "grace/grace_antarctica_mass.txt"), os.path.join(SRC, "grace/grace_greenland_mass.txt"),
           os.path.join(SRC, "glambie/calendar_years/0_global.csv"),
           os.path.join(SRC, "glambie/calendar_years/5_greenland_periphery.csv"),
           os.path.join(SRC, "noaa/noaa_thermosteric_w0-2000m_yearly.dat"),
           os.path.join(SRC, "noaa/nasa_gmsl_annual.csv"), LWS_GRACE_CSV]
OUT_PROV = os.path.join(_ARGS.work, "outputs/recalib_targets_ext_provenance.txt")
with open(OUT_PROV, "w") as _f:
    _f.write(
        f"built            {dt.datetime.now():%Y-%m-%d %H:%M:%S}\n"
        f"builder          tools/targets/prep_recalib_targets_ext.py (Ladrillo commit {_ARGS.commit})\n"
        f"file             {os.path.basename(OUT_CSV)}\n"
        f"sha256           {_sha(OUT_CSV)}\n"
        f"ais_source       {AIS_SOURCE}\n"
        f"reproduces       {AIS_SOURCE_REPRODUCES[AIS_SOURCE]}\n"
        f"ais_composition  Frederikse {FIT_Y0}-2018; GRACE {SPLICE_FROM['ais']}-{EXT_Y1} offset-matched over {OVERLAP['ais']}\n"
        f"span             {FIT_Y0}-{EXT_Y1}\n"
        f"baseline         {BASE_Y0}-{BASE_Y1}\n"
        f"lws_grace_tag    {LWS_GRACE_TAG}\n"
        f"units            cm SLE\n"
        f"randomness       none (no RNG is used)\n"
        + "".join(f"input            {_sha(p)}  {os.path.relpath(p, os.path.dirname(SRC)) if p.startswith(SRC) else os.path.relpath(p, _ARGS.work)}\n" for p in _inputs))
print(f"Wrote {OUT_PROV}")
