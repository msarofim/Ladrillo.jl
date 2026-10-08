"""
build_emissions_v160_cmip7harm_vv.py

van Vuuren ScenarioMIP-CMIP7 marker emissions on the **CMIP7 history** basis, so
the van Vuuren arm runs calib 1.6.0 with the inputs 1.6.0 was conditioned on.
One file per marker, same construction and same `cmip7harm` provenance infix as
the SSP builders.

WHY THIS EXISTS -- THE DEFECT IT FIXES
--------------------------------------
Every v160 van Vuuren driver -- `build_fair_cube_vv_v160.py`,
`build_fair_pulse_vv_v160.py`, `run_fair_vv_for_ccx.py`,
`run_fair_vanvuuren_ext_erf.py --calib-dir calibration_v160` -- reads
`data/vanvuuren/spliced_ext_{harmonized,zenodo}/`, whose 1750.5-2023.5 block is
**SMITH 2024 history**, because those files were built for calib 1.4.5 and were
never rebased. So the 1.6.0 posterior was being paired with v1.4.5 history: the
rung-R1 input mismatch the migration ladder exists to measure
(`scripts/diag_calib160_igcc_ladder.py`, memory `calib160_rffsp_ladder`).

MEASURED on marker M and marker H, 841 configs, holding the posterior, the
future and the marker's own prescribed forcing fixed and moving ONLY the history
block (2026-10-02, `notes/handoff_2026-10-02_calib160_future_emissions.md`):

    model 1850-1900 baseline GMST    -0.0619 K  ->  -0.1269 K  (0.065 K COOLER)
    GMST 2015-2024 bias vs IGCC      -0.1092 K  ->  -0.0608 K
    RMSE 1850-2023 vs IGCC            0.1317 K  ->   0.1217 K
    GMST 2100, anomaly                 +0.054 K   (M and H alike)
    GMST 2100, ABSOLUTE level          -0.010 K
    total ERF 2024 / 2100 / 2150      -0.022 / -0.005 / +0.005 W/m2
    P(ERF>8.5) @2150, marker H        42.212 %  ->  42.212 %  (UNCHANGED)

⛔ READ THE SECOND AND SIXTH ROWS TOGETHER. The basis does not make the model
warmer; it makes the 1850-1900 REFERENCE 0.065 K cooler, which is what closes
the gap to IGCC -- and IGCC is itself a 1850-1900 anomaly, so that is the
quantity fair-calibrate constrains and the improvement is real. The absolute
warming level at 2100 moves -0.010 K and absolute forcing moves <= 0.022 W/m2,
so NO absolute-threshold statistic changes.
⛔ A FIRST PASS CLAIMED 'total ERF 2100 +0.116 W/m2' AND 'P(ERF>8.5) @2150
45.42 -> 51.84 %'. BOTH ARE RETRACTED (same day). They came from re-referencing
total ERF to each arm's own 1850-1900 mean before comparing it to 8.5 W/m2.
FaIR already reports forcing relative to 1750 -- median total ERF at 1750 is
exactly 0.00000 on both bases -- so that double-references it, and the CMIP7
basis's 1850-1900 mean ERF is 0.098 W/m2 lower, which the second referencing
injected as a shift that is not in the forcing.

THE CONSTRUCTION, AND WHY IT IS THE SSP CONVENTION UNCHANGED
------------------------------------------------------------
CMIP7 historical 1750.5-2023.5 **verbatim** from the 1.6.0 Zenodo record, then
the marker's own future from 2024.5, re-harmonized at 2023.5 by a PER-SPECIES
multiplicative scale = hist(2023.5) / vv(2023.5). Identical in form to
`build_emissions_v160_ssp_rcmip.py`, and gated the same way: the scenario's own
year-on-year step across the join is preserved to float64 (`MAX_JOIN_STEP_RESIDUAL`).

The input is the ALREADY-SPLICED marker file, not the raw IIASA/Zenodo sources,
and this script adds exactly one join, at the front. The front scale is
non-negative for every species, so it is sign-safe on CO2 FFI and CO2 AFOLU where
they cross zero after 2100 (zenodo: x1.0014 and x1.0000).

⛔ THE DEFAULT VARIANT IS `zenodo` SINCE 2026-10-08 (Marcus: "go with b").
Until then it was `harmonized`: the IIASA prerelease to 2100, then the published
extension ADDED to the prerelease 2100 level for CO2
(`ours(2100) + [ext(t) - ext(2100)]`). That is exact only where the prerelease
and the published Zenodo v1.1.1 file agree at 2100, which was verified for H
only. Medium-to-Low was revised before publication: its CO2 at 2100 is -9.23
GtCO2/yr in the prerelease and +1.26 published. So the `harmonized` ML tail sat
10.49 GtCO2/yr low for 200 years: -2,098 GtCO2 of removal not in the scenario,
CO2 below pre-industrial by 2300, and 0.73 K too cold at 2300 (L -167 and LN
-112 GtCO2; H, HL, M and VL within 27). Measured in SLR-RFF-BRICK
`python/diag_vv_zenodo_ab/` (CHANGELOG 2026-10-08j). The continuity gate could
not see it: continuity at 2100 held BY CONSTRUCTION.
`zenodo` takes the WHOLE scenario period 2024.5-2300.5 from the published file,
so there is no 2100 join at all. Its history vintage differs from the 1.6.0
record (ozone-depleting gases 5-25% lower), so the front scale is no longer the
identity; that choice (b: keep the 1.6.0 history, scale the future) was measured
against (a: take Zenodo's own history) on ML at <= 0.009 K. See
SCALE_DEVIATION_OVERRIDE for the three gases whose scale leaves the 25% band.

THE ALTERNATIVE WE DID NOT TAKE, AND WHY IT IS STILL USEFUL
-----------------------------------------------------------
`data/vanvuuren/spliced_ext_fullnative/` already carries a CMIP7-ERA history
(Nicholls et al. 2025) over the whole 1750-2300 span with NO join at all. It is
not this script's basis for one reason: its history is *near* the 1.6.0 record's
but not the record's own bytes (CH4 2000.5: 307.4514 vs 307.4902 Mt/yr; CO2 FFI
2015.5: 35.9590 vs 36.0139 Gt/yr), and the standing convention since 2026-08-28
is CMIP7 history **verbatim**, which is what makes the splice auditable against
the Zenodo file. `fullnative` remains the right arm for "is the gap emissions or
calibration vintage", and a `fullnative`-vs-`cmip7harm` difference is a clean
measure of what the history VINTAGE within CMIP7 is worth.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_emissions_v160_cmip7harm_vv.py @ 8984a13. The construction is
copied UNCHANGED; only the I/O plumbing differs:
  * --hist-file, --spliced-ext-root and --out-dir are REQUIRED. The van Vuuren
    template is read from <spliced-ext-root>/spliced_ext_<variant>/ (was
    FaIRtoFrEDI/data/vanvuuren/), history from --hist-file (was v160_cache/),
    and the output goes to --out-dir (was calibration_v160_prod/).
Normally driven by build_forcing.py.
"""

from pathlib import Path
import argparse
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fair_units import fair_canonical_unit, report as unit_report

# -- Named constants (every label and filename below derives from these) ------
CALIB_VERSION = "1.6.0"
HIST_LABEL = "CMIP7 historical 1750-2023 (Zenodo 18828694)"
BASIS_SLUG = "cmip7harm"            # the provenance infix; see module docstring
SCEN_PREFIX = "vv"                  # output key = SCEN_PREFIX + MARKER, e.g. vvM

HIST_FILE = None             # set from --hist-file in main()
SPLICED_EXT_ROOT = None      # set from --spliced-ext-root in main()

# ── THE SCENARIO VARIANT, AND HOW IT REACHES THE FILENAME ─────────────────────
# `harmonized` and `native` share the IIASA prerelease 2024-2100 and differ only
# after 2100.5; `zenodo` is the published file over the whole 2024.5-2300.5. For
# the prerelease variants the front rescale is the identity to ~1e-8; for
# `zenodo` it is not (see SCALE_DEVIATION_OVERRIDE). The history basis and the
# variant are independent axes. Both must be readable from the filename, or two
# different objects collide at one path.
#   DEFAULT_VARIANT  -> emissions_v160_cmip7harm_vv<MARKER>.csv
#   any other        -> emissions_v160_cmip7harm_<variant>_vv<MARKER>.csv
# The unsuffixed spelling meaning the DEFAULT is a stated convention, the same
# one calibration_v160_prod already uses for volcanic_solar.csv (unsuffixed =
# marker-free). It is stated in both READMEs and printed on every run.
# ⚠ Variants COEXIST SAFELY in one directory, verified rather than assumed: the
# cube and pulse drivers glob `emissions_v160_cmip7harm_vv*.csv`, and `vv` has to
# sit immediately after `cmip7harm_`, so `..._cmip7harm_zenodo_vvM.csv` does NOT
# match (checked with fnmatch, 2026-10-02). An earlier version of this comment
# claimed the opposite -- that a second variant would make resolve_marker() match
# two files per marker and exit. It would not; it is simply ignored.
DEFAULT_VARIANT = "zenodo"          # the cube drivers' basis since 2026-10-08 (was "harmonized"; see docstring)
VARIANTS = ("harmonized", "zenodo", "native")

HIST_LAST_YEAR_INT = 2023           # CMIP7 history ends here = the join anchor
DROP_SPECIES = {"Halon-1202"}       # in CMIP7, not a species in 1.4.5 or 1.6.0

# A harmonized series must preserve the scenario's OWN year-on-year step across
# the join; anything larger is a broken anchor, not a property of the data.
MAX_JOIN_STEP_RESIDUAL = 1e-9
# The front rescale is a level correction between two historical reconstructions
# of the same world. A species moving more than this is a unit or definition
# mismatch, not a reconstruction difference -- stop rather than scale it.
MAX_ABS_SCALE_DEVIATION = 0.25
# MEASURED exceptions, not loosened tolerances (2026-10-08; memory
# `zenodo_vv_history_vintage`). The Zenodo v1.1.1 file carries its OWN history
# vintage, and for these gases the two reconstructions genuinely differ by more
# than 25% at 2023 -- not a unit error (a unit error is x1000):
#   CCl4        x1.257  Zenodo's history runs 10-22% below the 1.6.0 record in
#                       every decade since 1960, so 2023 is not an outlier
#   CFC-115     x0.181  near zero and noisy; BOTH records hold 0 in 2021-2022
#                       and a non-zero 2023 (a data oddity in both)
#   Halon-2402  x0      the 1.6.0 record has 0 in 2023, so the scale zeroes the
#                       future (published 2023: 0.28 kt/yr)
# Each value is the largest |scale - 1| the species may take; it applies to every
# variant, and every other species keeps MAX_ABS_SCALE_DEVIATION.
SCALE_DEVIATION_OVERRIDE = {"CCl4": 0.30, "CFC-115": 0.85, "Halon-2402": 1.0}

# Declared, not parsed out of the prose scenario name: "High-to-Low - SSP5
# (Marker)" -> HL is a READING, and a reading belongs in a constant a human can
# check. Same table as build_fair_cube_vv_v160.py, by design.
SCENARIO_TO_MARKER = {
    "High - SSP3 (Marker)":            "H",
    "High-to-Low - SSP5 (Marker)":     "HL",
    "Medium - SSP2 (Marker)":          "M",
    "Medium-to-Low - SSP2 (Marker)":   "ML",
    "Low - SSP2 (Marker)":             "L",
    "Low-to-Negative - SSP2 (Marker)": "LN",
    "Very Low - SSP1 (Marker)":        "VL",
}


def emis_dir(variant):
    return SPLICED_EXT_ROOT / f"spliced_ext_{variant}"


def marker_files(variant):
    """marker code -> its spliced emissions CSV, gating the whole table on disk.

    The mapping must be a bijection over exactly the seven markers and every
    scenario found on disk must be in the table, so an unrecognised scenario is
    an error at startup rather than a cube built on the wrong forcing.
    """
    d = emis_dir(variant)
    if len(set(SCENARIO_TO_MARKER.values())) != len(SCENARIO_TO_MARKER):
        sys.exit("ERROR: SCENARIO_TO_MARKER is not a bijection")
    csvs = sorted(d.glob("*.csv"))
    if not csvs:
        sys.exit(f"ERROR: no emissions CSVs in {d}")
    out = {}
    for p in csvs:
        scen = str(pd.read_csv(p, nrows=1)["scenario"].iloc[0])
        if scen not in SCENARIO_TO_MARKER:
            sys.exit(f"ERROR: {p.name} holds unknown scenario {scen!r}; add it to "
                     f"SCENARIO_TO_MARKER or remove the file")
        m = SCENARIO_TO_MARKER[scen]
        if m in out:
            sys.exit(f"ERROR: marker {m} matched both {out[m][0].name} and {p.name}")
        out[m] = (p, scen)
    missing = sorted(set(SCENARIO_TO_MARKER.values()) - set(out))
    if missing:
        sys.exit(f"ERROR: {d.name} has no file for marker(s) {missing}")
    return out


def load_history():
    dfh = pd.read_csv(HIST_FILE)
    dfh = dfh[~dfh["variable"].isin(DROP_SPECIES)].set_index("variable")
    ycols = [c for c in dfh.columns if str(c).strip().isdigit()]
    if max(int(c) for c in ycols) != HIST_LAST_YEAR_INT:
        sys.exit(f"ERROR: {HIST_FILE.name} does not end at {HIST_LAST_YEAR_INT}")
    return dfh, ycols, [f"{int(c)}.5" for c in ycols]


def build_one(marker, path, scen, dfh, ycols, hist_cols):
    dv = pd.read_csv(path).set_index("variable")
    vcols = [c for c in dv.columns if str(c).replace(".", "").isdigit()]
    a_col = f"{HIST_LAST_YEAR_INT}.5"
    if a_col not in vcols:
        sys.exit(f"ERROR: {path.name} has no {a_col} column to anchor on")
    fut_cols = [c for c in vcols if float(c) > HIST_LAST_YEAR_INT + 0.5]

    missing = sorted(set(dv.index) - set(dfh.index))
    if missing:
        sys.exit(f"ERROR: CMIP7 history lacks {path.name} species {missing}")

    rows, renamed, scales = [], {}, {}
    worst, worst_sp = 0.0, None
    for sp in dv.index:
        h = dfh.loc[sp, ycols].to_numpy(float)
        h_a, v_a = float(h[-1]), float(dv.loc[sp, a_col])
        scale = (h_a / v_a) if abs(v_a) > 1e-12 else (1.0 if abs(h_a) < 1e-12 else 0.0)
        band = SCALE_DEVIATION_OVERRIDE.get(sp, MAX_ABS_SCALE_DEVIATION)
        if abs(scale - 1.0) > band:
            sys.exit(f"ERROR: {sp} needs a x{scale:.4f} front rescale, beyond the "
                     f"{band:.0%} band two reconstructions of the "
                     f"same history should differ by. Suspect units or a definition "
                     f"change, not a reconstruction difference.")
        scales[sp] = scale
        fut = dv.loc[sp, fut_cols].to_numpy(float) * scale

        src = str(dv.loc[sp, "unit"]).strip()
        unit = fair_canonical_unit(sp, src)
        if unit != src:
            renamed[sp] = (src, unit)

        # The source's OWN prose scenario name, verbatim -- "Medium - SSP2
        # (Marker)", not "vvM". The existing VV drivers identify a marker by
        # reading this column and looking it up in SCENARIO_TO_MARKER, so
        # renaming it here would make every one of them reject the file. The
        # BASIS belongs in the filename (BASIS_SLUG); the SCENARIO belongs here.
        row = {"scenario": scen, "region": "World", "variable": sp, "unit": unit}
        row.update(dict(zip(hist_cols, h)))
        row.update(dict(zip(fut_cols, fut)))
        rows.append(row)

        v_b = float(dv.loc[sp, fut_cols[0]])
        if abs(h_a) > 1e-12 and abs(v_a) > 1e-12:
            d = abs(float(fut[0]) / h_a - v_b / v_a)
            if d > worst:
                worst, worst_sp = d, sp

    out = pd.DataFrame(rows, columns=["scenario", "region", "variable", "unit",
                                      *hist_cols, *fut_cols])
    yrs = [float(c) for c in [*hist_cols, *fut_cols]]
    if not all(b > a for a, b in zip(yrs, yrs[1:])):
        sys.exit(f"ERROR: {marker} year columns not strictly increasing -- FaIR "
                 f"will raise NonMonotonicError")
    if out.isna().any().any():
        sys.exit(f"ERROR: NaN in {marker} output")
    if worst >= MAX_JOIN_STEP_RESIDUAL:
        sys.exit(f"ERROR: {marker} harmonization broke the scenario's own step at "
                 f"the join ({worst:.3e}, {worst_sp})")

    # The history block must be the record's own bytes, not a near copy: that is
    # the whole claim the `cmip7harm` infix makes. Compare what we WROTE against
    # what we READ, so an edit to the assembly above cannot pass this silently.
    reread = out.set_index("variable")[hist_cols].to_numpy(float)
    recorded = dfh.loc[out["variable"], ycols].to_numpy(float)
    if not np.array_equal(reread, recorded):
        sys.exit(f"ERROR: {marker} history block is not identical to "
                 f"{HIST_FILE.name}; the cmip7harm label would be false")

    # NO '#' provenance header: fair.io.fill_from_csv does not skip comments, so
    # one would make the file unreadable by the driver that consumes it.
    # Provenance lives in the BASIS_SLUG infix and in calibration_v160_prod/README.md.
    return out, scales, worst, worst_sp, renamed, scen, len(fut_cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", nargs="*", default=None,
                    help=f"subset of {sorted(set(SCENARIO_TO_MARKER.values()))}; "
                         f"default all seven")
    ap.add_argument("--variant", default=DEFAULT_VARIANT, choices=VARIANTS,
                    help=f"which spliced_ext_* tail convention to rebase "
                         f"(default {DEFAULT_VARIANT})")
    ap.add_argument("--hist-file", type=Path, required=True,
                    help="historical_emissions_1750-2023_cmip7.csv (Zenodo 18828694)")
    ap.add_argument("--spliced-ext-root", type=Path, required=True,
                    help="directory holding spliced_ext_<variant>/ (written by "
                         "build_vanvuuren_extended_splice.py)")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    global HIST_FILE, SPLICED_EXT_ROOT
    HIST_FILE, SPLICED_EXT_ROOT = args.hist_file, args.spliced_ext_root
    args.out_dir.mkdir(parents=True, exist_ok=True)

    files = marker_files(args.variant)
    want = args.markers or sorted(files)
    unknown = sorted(set(want) - set(files))
    if unknown:
        sys.exit(f"ERROR: unknown marker(s) {unknown}; have {sorted(files)}")

    dfh, ycols, hist_cols = load_history()
    print(f"-- van Vuuren emissions on the {BASIS_SLUG} basis, calib {CALIB_VERSION} --")
    print(f"  history : {HIST_LABEL}")
    print(f"  future  : data/vanvuuren/spliced_ext_{args.variant}/ from "
          f"{HIST_LAST_YEAR_INT + 1}.5")
    print(f"  join    : {HIST_LAST_YEAR_INT}.5, per-species multiplicative")
    print(f"  markers : {' '.join(want)}")
    print(f"  tail    : spliced_ext_{args.variant}"
          + ("  (the DEFAULT; written without a variant infix)"
             if args.variant == DEFAULT_VARIANT
             else f"  (NON-default; filename carries '{args.variant}_')"))

    for m in want:
        path, scen = files[m]
        out, scales, worst, worst_sp, renamed, scen, n_fut = build_one(
            m, path, scen, dfh, ycols, hist_cols)
        unit_report(renamed)
        vtag = "" if args.variant == DEFAULT_VARIANT else f"{args.variant}_"
        dest = (args.out_dir
                / f"emissions_v160_{BASIS_SLUG}_{vtag}{SCEN_PREFIX}{m}.csv")
        out.to_csv(dest, index=False)
        # The variant travels WITH the file (a CSV comment would break fill_from_csv), so the cube
        # driver's provenance is read off what was built, not asserted from today's default.
        dest.with_suffix(".variant").write_text(args.variant + "\n")
        s = pd.Series(scales)
        print(f"\n  {m} ({scen}): {out.shape[0]} species x "
              f"{len(hist_cols) + n_fut} timepoints "
              f"({hist_cols[0]} .. {out.columns[-1]})")
        print(f"    history block identical to {HIST_FILE.name}  PASS")
        print(f"    join gate worst step residual {worst:.2e} ({worst_sp})  PASS")
        # ⚠ PRINTED IN SCIENTIFIC NOTATION DELIBERATELY. At 4 decimal places every
        # scale reads "x1.0000" and "max |dev| 0.0000", which invites the reader to
        # conclude the rescale is EXACTLY the identity. It is not: the van Vuuren
        # future was itself harmonized onto a CMIP7-era history at its own 2023
        # join, so hist(2023.5)/vv(2023.5) is 1 to within ~1e-8 but not 1 -- and
        # "suspiciously uniform" is a bug signal that should not be manufactured
        # by a format string.
        print(f"    front rescale: median x{s.median():.10f}, max |dev| "
              f"{float(s.sub(1).abs().max()):.3e} ({s.sub(1).abs().idxmax()})")
        print(f"    wrote {dest.name}")

    print(f"\nNOTE: the '{BASIS_SLUG}' infix is the provenance marker. The van "
          f"Vuuren files under data/vanvuuren/spliced_ext_*/ are the SMITH 2024 "
          f"history basis and belong to the calib 1.4.5 lineage; pair these with "
          f"calib {CALIB_VERSION}, and pair each marker with its OWN "
          f"volcanic_solar_<MARKER>.csv.")


if __name__ == "__main__":
    main()
