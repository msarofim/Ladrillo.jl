"""
build_emissions_v160_ssp245.py

SSP2-4.5 emissions for the fair-calibrate **1.6.0** production stack: CMIP7
historical 1750-2023 spliced to chrisroadmap SSP2-4.5, harmonized at the join.

This is the 1.6.0 analogue of `build_emissions_v145_ssp245.py` and follows the
SAME splice convention Marcus set on 2026-05-22 -- per-species MULTIPLICATIVE
scaling of the future so the series is continuous at the join -- with two
changes decided 2026-08-28:

  1. HISTORY is CMIP7 `historical_emissions_1750-2023_cmip7.csv` (the file
     fair-calibrate 1.6.0 was conditioned on) instead of Smith 2024.
  2. The JOIN moves 2020.5 -> 2023.5, taking all three extra years of observed
     history. `diag_calib160_igcc_ladder.py` measured this component at
     +0.05 K of historical GMST fit, which is why it is not optional.

Post-2100 extension is UNCHANGED from the v145 convention: per-species r(t)
lifted from `AR6_emissions_ssp245_1750_2300.csv`, with FossilCO2's r(t) as the
Group-A template for species AR6 does not carry or whose AR6 2100 is ~0.

WHY THE SPECIES MAPS ARE PARSED, NOT COPIED: `AR6_TO_V145`,
`V145_ORPHAN_SPECIES` and `GROUP_A_TEMPLATE_AR6` are read out of
`build_emissions_v145_ssp245.py` with `ast`, so the two builders cannot drift
apart. The v145 script is frozen provenance and is never imported or executed
(it is a flat script; importing it would run it).

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_emissions_v160_ssp245.py @ 8984a13. The construction is copied
UNCHANGED; only the I/O plumbing differs:
  * the three inputs are REQUIRED flags (--hist-file, --chris-file, --ar6-file)
    and the output goes to --out-dir (was FaIRtoFrEDI/calibration_v160).
  * the species maps are no longer parsed out of the PRIVATE research-repo
    script build_emissions_v145_ssp245.py; they are vendored verbatim in
    v145_species_maps.py. Pass --v145-builder PATH to re-parse that script and
    REFUSE to run if the vendored copy has drifted from it.
  * the closing "CMIP7-history vs the v145 file it replaces" printout read a
    research-repo file; it now runs only with --compare-v145-file. It is a
    printout after the write and never touched the output.
Normally driven by build_forcing.py.
"""

from pathlib import Path
import argparse
import ast
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
# The two-warrant check that LIVED HERE is now shared, so the four v160
# builders cannot drift apart on it. UNIT_ALIASES moved with it.
from fair_units import UNIT_ALIASES, fair_canonical_unit

# ── Named constants: every label and filename below derives from these ───────
CALIB_VERSION = "1.6.0"
HIST_SOURCE_LABEL = "CMIP7 historical 1750-2023 (Zenodo 18828694)"
FUTURE_SOURCE_LABEL = "chrisroadmap MESSAGE-GLOBIOM SSP2-4.5"

# Input paths come from the CLI (see main()); the output NAME is fixed here.
HIST_FILE = CHRIS_FILE = AR6_FILE = V145_BUILDER = None
OUT_NAME = "emissions_v160_cmip7_ssp245_harmonized.csv"
OUT_HARMON = None

HIST_LAST_YEAR_INT = 2023     # CMIP7 historical ends here; the harmonization anchor
CHRIS_LAST_YEAR_INT = 2100    # chrisroadmap MESSAGE-GLOBIOM file ends here
END_YEAR_INT = 2300           # final timepoint 2300.5

CHRIS_SCENARIO = "MESSAGE-GLOBIOM___SSP2-45"   # plain SSP2-4.5, not the CH4-pledge variants

# Present in BOTH the CMIP7 historical file and chrisroadmap, but NOT a species
# in fair-calibrate 1.4.5 or 1.6.0 (1.6.0's species list is 1.4.5's plus
# Irrigation). Dropped exactly as build_emissions_v145_ssp245.py drops it, so
# the output's species set matches what FaIR will define.
DROP_SPECIES = {"Halon-1202"}

# A harmonization that is doing its job leaves NO jump at the join. Anything
# above this is a bug in the anchor, not a property of the data.
MAX_JOIN_JUMP_PCT = 1e-6

# CMIP7 and chrisroadmap spell two units differently for the SAME quantity. The
# unit gate below stays STRICT -- a real unit mismatch would silently rescale the
# multiplicative anchor -- so equivalences are enumerated here rather than the
# check being relaxed. Both are verified label-only:
#   NOx          CMIP7 "Mt NO2/yr" == chris "Mt NOx/yr". Smith 2024 also says
#                "Mt NO2/yr" and build_emissions_v145_ssp245.py already splices
#                the two, and 2023 values agree in magnitude (127.5 vs 146.3),
#                which MtN vs MtNO2 (a factor 46/14) could not.
#   HFC-4310mee  CMIP7 "kt HFC4310/yr" == chris "kt HFC4310mee/yr" -- the same
#                species with the isomer suffix dropped.


# ⛔⛔ THE ALIAS SET ABOVE ONLY MAKES THE *GATE* PERMISSIVE. Until 2026-09-24 it
# did NOT normalise what is WRITTEN: the row below took its unit straight from
# the CMIP7 history frame, so HFC-4310mee shipped as `kt HFC4310/yr` -- a
# spelling FaIR does not know. FaIR then takes its NOVEL-COMPOUND branch
# (fair/io/fill_from.py:95-101) and writes its MODULE-LEVEL unit tables:
# `desired_emissions_units["HFC-4310mee"] = "kt HFC4310/yr"` and
# `compound_convert["HFC4310"] = {...}`. The NEXT scenario's fill_from_rcmip()
# in the same process then dies with `KeyError: 'HFC4310'`.
# ⚠⚠ THE TELL IS THAT SCENARIO 1 PASSES AND SCENARIO 2 FAILS, so a
# single-scenario smoke test cannot see it -- which is how it shipped.
#
# ⛔ THE CANONICAL SPELLING IS FaIR's OWN, AND IT IS NOT ALWAYS THE ALIAS PAIR'S
# SECOND ELEMENT. Measured 2026-09-24: FaIR wants `kt HFC4310mee/yr` (the second)
# but `Mt NO2/yr` for NOx (the FIRST -- the CMIP7 spelling). A blind normalise to
# one side of the pair would have fixed one species and broken the other, so the
# spelling is READ OUT OF FaIR rather than hand-maintained here.
# ⚠ Relabelling without rescaling is only sound at a conversion factor of
# exactly 1.0. Asserted, not assumed -- both are 1.0 in FaIR's own table.


def load_v145_maps():
    """The species maps, vendored verbatim; cross-checked against the v145
    builder's own definitions when --v145-builder is given."""
    import v145_species_maps as vm
    vendored = (vm.AR6_TO_V145, vm.V145_ORPHAN_SPECIES, vm.GROUP_A_TEMPLATE_AR6)
    if V145_BUILDER is None:
        print("  maps    : vendored v145_species_maps.py (no --v145-builder "
              "cross-check)")
        return vendored
    parsed = _parse_v145_maps()
    if parsed != vendored:
        sys.exit(f"ERROR: v145_species_maps.py differs from {V145_BUILDER.name}"
                 f" -- the vendored copy has drifted; fix before regenerating.")
    print(f"  maps    : vendored copy == {V145_BUILDER.name}  PASS")
    return vendored


def _parse_v145_maps():
    """Parse the species maps out of the v145 builder without executing it."""
    tree = ast.parse(V145_BUILDER.read_text())
    want = {"AR6_TO_V145", "V145_ORPHAN_SPECIES", "GROUP_A_TEMPLATE_AR6"}
    found = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            t = node.targets[0]
            if isinstance(t, ast.Name) and t.id in want:
                found[t.id] = ast.literal_eval(node.value)
    missing = want - set(found)
    if missing:
        sys.exit(f"ERROR: {V145_BUILDER.name} no longer defines {sorted(missing)} "
                 f"-- the two builders have diverged; fix before regenerating.")
    return found["AR6_TO_V145"], found["V145_ORPHAN_SPECIES"], found["GROUP_A_TEMPLATE_AR6"]


def main():
    global HIST_FILE, CHRIS_FILE, AR6_FILE, V145_BUILDER, OUT_HARMON
    ap = argparse.ArgumentParser()
    ap.add_argument("--hist-file", type=Path, required=True,
                    help="historical_emissions_1750-2023_cmip7.csv (Zenodo 18828694)")
    ap.add_argument("--chris-file", type=Path, required=True,
                    help="chrisroadmap CH4Pledge-rcmip-emissions-ssp245-CLE-MFR-filled.csv")
    ap.add_argument("--ar6-file", type=Path, required=True,
                    help="AR6_emissions_ssp245_1750_2300.csv (post-2100 r(t) shape)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--v145-builder", type=Path, default=None,
                    help="research-repo build_emissions_v145_ssp245.py, to "
                         "cross-check the vendored species maps")
    ap.add_argument("--compare-v145-file", type=Path, default=None,
                    help="emissions_v145_ssp245_harmonized.csv, for the closing "
                         "printout only")
    args = ap.parse_args()
    HIST_FILE, CHRIS_FILE, AR6_FILE = args.hist_file, args.chris_file, args.ar6_file
    V145_BUILDER = args.v145_builder
    OUT_HARMON = args.out_dir / OUT_NAME

    AR6_TO_V145, ORPHAN_SPECIES, GROUP_A_TEMPLATE_AR6 = load_v145_maps()
    print("── SSP2-4.5 emissions for calib 1.6.0 ───────────────────────────────")
    print(f"  history : {HIST_SOURCE_LABEL}")
    print(f"  future  : {FUTURE_SOURCE_LABEL}")
    print(f"  join    : {HIST_LAST_YEAR_INT}.5, per-species multiplicative")
    print(f"  maps    : {len(AR6_TO_V145)} AR6 names, {len(ORPHAN_SPECIES)} orphans")

    # ── history: CMIP7, integer year columns -> FaIR timepoints ──────────────
    dfh = pd.read_csv(HIST_FILE)
    ycols = [c for c in dfh.columns if str(c).strip().isdigit()]
    hist_years = [int(c) for c in ycols]
    if max(hist_years) != HIST_LAST_YEAR_INT:
        sys.exit(f"ERROR: {HIST_FILE.name} ends {max(hist_years)}, expected "
                 f"{HIST_LAST_YEAR_INT}")
    hist_cols = [f"{y}.5" for y in hist_years]
    hist = pd.DataFrame(dfh[ycols].to_numpy(float), columns=hist_cols)
    hist.insert(0, "variable", dfh["variable"].to_numpy())
    hist.insert(1, "unit", dfh["unit"].to_numpy())
    hist = hist[~hist["variable"].isin(DROP_SPECIES)].reset_index(drop=True)
    print(f"\n  history: {len(hist)} species, {hist_cols[0]} .. {hist_cols[-1]}")

    # ── future: chrisroadmap ─────────────────────────────────────────────────
    dfc = pd.read_csv(CHRIS_FILE)
    dfc = dfc[dfc["scenario"] == CHRIS_SCENARIO].copy()
    dfc = dfc[~dfc["variable"].isin(DROP_SPECIES)].copy()
    meta = ("scenario", "region", "variable", "unit")
    chris_future_cols = [c for c in dfc.columns if c not in meta
                         and HIST_LAST_YEAR_INT + 0.5 < float(c) <= CHRIS_LAST_YEAR_INT + 0.5]
    print(f"  future : {len(dfc)} species, {chris_future_cols[0]} .. "
          f"{chris_future_cols[-1]}")

    hist_sp, chris_sp = set(hist["variable"]), set(dfc["variable"])
    if hist_sp != chris_sp:
        sys.exit(f"ERROR: species mismatch. CMIP7 only: {sorted(hist_sp - chris_sp)}; "
                 f"chris only: {sorted(chris_sp - hist_sp)}")
    # Units must match species by species or the multiplicative anchor silently
    # rescales between two different unit systems.
    hu = hist.set_index("variable")["unit"]
    cu = dfc.set_index("variable")["unit"]
    def units_ok(x, y):
        x, y = str(x).strip(), str(y).strip()
        return x == y or (x, y) in UNIT_ALIASES or (y, x) in UNIT_ALIASES
    bad = [v for v in sorted(hist_sp) if not units_ok(hu[v], cu[v])]
    aliased = [v for v in sorted(hist_sp)
               if str(hu[v]).strip() != str(cu[v]).strip() and units_ok(hu[v], cu[v])]
    if bad:
        sys.exit(f"ERROR: unit mismatch between history and future for {bad}")
    print(f"  species: {len(hist_sp)} matched, units agree"
          + (f" ({len(aliased)} via known label alias: {aliased})" if aliased else ""))

    # ── post-2100 r(t) from the AR6 wide file (v145 convention, unchanged) ────
    ar6 = pd.read_csv(AR6_FILE).set_index("Year")
    ext_years = list(range(CHRIS_LAST_YEAR_INT + 1, END_YEAR_INT + 1))
    ar6_2100 = ar6.loc[CHRIS_LAST_YEAR_INT]
    template_r = (ar6.loc[ext_years, GROUP_A_TEMPLATE_AR6].to_numpy(float)
                  / float(ar6_2100[GROUP_A_TEMPLATE_AR6]))
    r_dict, degenerate = {}, []
    for ar6_name, v_name in AR6_TO_V145.items():
        if v_name is None:
            continue
        denom = float(ar6_2100[ar6_name])
        if abs(denom) < 1e-12:
            r_dict[v_name] = template_r.copy()
            degenerate.append(v_name)
        else:
            r_dict[v_name] = ar6.loc[ext_years, ar6_name].to_numpy(float) / denom
    for sp in ORPHAN_SPECIES:
        r_dict[sp] = template_r.copy()
    missing = hist_sp - set(r_dict)
    if missing:
        sys.exit(f"ERROR: no post-2100 r(t) for {sorted(missing)}")
    print(f"  r(t)   : {len(ext_years)} years; {len(degenerate)} degenerate + "
          f"{len(ORPHAN_SPECIES)} orphan species use the "
          f"{GROUP_A_TEMPLATE_AR6} template")

    # ── build ────────────────────────────────────────────────────────────────
    ext_cols = [f"{y}.5" for y in ext_years]
    anchor = f"{HIST_LAST_YEAR_INT}.5"
    hist_i = hist.set_index("variable")
    chris_i = dfc.set_index("variable")

    rows, scales = [], {}
    for v in hist["variable"]:
        h, c = hist_i.loc[v], chris_i.loc[v]
        h_a, c_a = float(h[anchor]), float(c[anchor])
        # Zero-denominator: a species the scenario retires. Keep it at zero
        # rather than dividing; carrying history's value forward would invent
        # emissions the scenario does not have.
        scale = (h_a / c_a) if abs(c_a) > 1e-12 else (1.0 if abs(h_a) < 1e-12 else 0.0)
        scales[v] = scale
        row = {"scenario": "all", "region": "World", "variable": v,
               "unit": fair_canonical_unit(v, h["unit"])}
        row.update({c_: float(h[c_]) for c_ in hist_cols})
        row.update({c_: float(c[c_]) * scale for c_ in chris_future_cols})
        c_2100 = float(c[f"{CHRIS_LAST_YEAR_INT}.5"])
        row.update({c_: c_2100 * scale * r_dict[v][i] for i, c_ in enumerate(ext_cols)})
        rows.append(row)

    all_cols = ["scenario", "region", "variable", "unit",
                *hist_cols, *chris_future_cols, *ext_cols]
    out = pd.DataFrame(rows, columns=all_cols)

    # ── gates ────────────────────────────────────────────────────────────────
    yrs = [float(c) for c in all_cols[4:]]
    if not all(b > a for a, b in zip(yrs, yrs[1:])):
        sys.exit("ERROR: year columns are not strictly increasing -- FaIR will "
                 "raise NonMonotonicError")
    # ⭐ THE GATE THAT WOULD HAVE CAUGHT THIS: every unit written must name a
    # compound FaIR ALREADY KNOWS, so the novel-compound branch is never taken.
    # Derived from the CONSUMER's own table rather than a list maintained here.
    from fair.structure.units import compound_convert as _cc
    unknown = sorted({(r.variable, r.unit) for r in out.itertuples()
                      if r.unit.split()[-1].split("/")[0] not in _cc})
    if unknown:
        sys.exit("ERROR: units name compounds FaIR does not know, so fill_from_csv "
                 "would mutate FaIR's MODULE-LEVEL unit tables and the NEXT "
                 f"scenario in the same process would fail: {unknown}")
    if out.isna().any().any():
        sys.exit("ERROR: NaN in output")

    i_a = all_cols.index(anchor)
    nxt = all_cols[i_a + 1]
    worst, worst_sp = 0.0, None
    for _, r in out.iterrows():
        a, b = float(r[anchor]), float(r[nxt])
        if abs(a) < 1e-12:
            continue
        # The harmonization makes the FUTURE meet HISTORY at the anchor; the
        # residual jump is the scenario's own year-on-year change, so compare
        # against the UNSCALED future's own step, not against zero.
        c_a, c_b = float(chris_i.loc[r["variable"]][anchor]), float(chris_i.loc[r["variable"]][nxt])
        if abs(c_a) < 1e-12:
            continue
        d = abs((b / a) - (c_b / c_a))
        if d > worst:
            worst, worst_sp = d, r["variable"]
    print(f"\n  join gate: worst |harmonized step - native step| = {worst:.3e} "
          f"({worst_sp})   {'PASS' if worst < MAX_JOIN_JUMP_PCT else 'FAIL'}")
    if worst >= MAX_JOIN_JUMP_PCT:
        sys.exit("ERROR: harmonization did not preserve the scenario's own step "
                 "at the join -- the anchor is wrong")

    OUT_HARMON.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_HARMON, index=False)
    print(f"\nWrote {OUT_HARMON}  ({out.shape})")

    print("\n=== how far the future moved (harmonization scale, "
          "history_2023.5 / chris_2023.5) ===")
    s = pd.Series(scales).sort_values()
    for v in list(s.index[:5]) + list(s.index[-5:]):
        print(f"  {v:>12s}  x{s[v]:.4f}")
    print(f"  median x{s.median():.4f}   "
          f"{(s.sub(1).abs() > 0.10).sum()} of {len(s)} species move >10%")

    if args.compare_v145_file is None:
        return
    print("\n=== key species, CMIP7-history vs the v145 file it replaces ===")
    old = pd.read_csv(args.compare_v145_file).set_index("variable")
    new = out.set_index("variable")
    for sp in ["CO2 FFI", "CO2 AFOLU", "CH4", "N2O", "Sulfur", "BC"]:
        line = f"  {sp:>10s} [{new.loc[sp,'unit']}]  "
        for y in ["2000.5", "2019.5", "2023.5", "2050.5", "2100.5"]:
            o = float(old.loc[sp, y]) if y in old.columns else float("nan")
            n = float(new.loc[sp, y])
            line += f"{y}: {o:.3g}->{n:.3g}   "
        print(line)


if __name__ == "__main__":
    main()
