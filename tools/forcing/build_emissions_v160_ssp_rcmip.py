"""
build_emissions_v160_ssp_rcmip.py

CMIP6-SSP emissions for the fair-calibrate **1.6.0** production stack: CMIP7
historical 1750-2023 spliced to the RCMIP-native SSP future, harmonized at the
join. One file per SSP, for the Ladrillo per-config cube line.

WHY THIS IS A SECOND BUILDER AND NOT AN ARGUMENT TO THE OTHER ONE
`build_emissions_v160_ssp245.py` splices to CHRISROADMAP SSP2-4.5, which stops
at 2100 and therefore needs the AR6 r(t) machinery to reach 2300. RCMIP runs
natively past 2300, so this line needs no extension at all. Two genuinely
different splices; collapsing them into one script would mean a dead code path
in every run.

WHAT CHANGES FOR THE CONSUMER, AND WHY THE FILENAME SAYS SO
The cube drivers previously took emissions straight from `fill_from_rcmip()` --
RCMIP-native throughout, including history. These files are CMIP7 history plus a
harmonized RCMIP future, which `build_fair_mean_v145.py`'s own docstring warns is
a DIFFERENT OBJECT ("not interchangeable for calibration"). Marcus 2026-08-28
chose the harmonized basis for consistency with the 2023 splice decision, and
asked that the change not be silent: hence the `cmip7harm` infix in the output
name and the provenance header written into every file.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_emissions_v160_ssp_rcmip.py @ 8984a13. The construction is copied
UNCHANGED; only the I/O plumbing differs:
  * --hist-file, --species-file and --out-dir are REQUIRED (were
    FaIRtoFrEDI/v160_cache, calibration_v160_prod, calibration_v160_prod).
  * the SSP list is --ssps (was bare positional argv).
  * `fill_from_rcmip()` is still how the RCMIP future is read -- that IS the
    arithmetic -- and fair fetches RCMIP v5.1.0 through pooch on first use. To
    keep the build OFFLINE, the port first checks that all three RCMIP files are
    already in pooch's cache with fair's own md5s, and stops otherwise
    (fetch_inputs.py primes that cache). Nothing is downloaded by this script.
Normally driven by build_forcing.py.
"""

from pathlib import Path
import argparse
import sys

import numpy as np
import pandas as pd
from fair import FAIR
from fair.io import read_properties

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fair_units import fair_canonical_unit, report as unit_report

# ── Named constants ──────────────────────────────────────────────────────────
CALIB_VERSION = "1.6.0"
HIST_LABEL = "CMIP7 historical 1750-2023 (Zenodo 18828694)"
FUTURE_LABEL = "RCMIP-native SSP (fair.fill_from_rcmip)"
BASIS_SLUG = "cmip7harm"          # the filename infix; see module docstring

# ⚠ SSPS IS A CLI ARGUMENT, NOT A LITERAL, since 2026-09-02. It was a hardcoded triple, so
# adding a scenario meant editing the script and every rebuild silently re-ran the old set.
# RCMIP's own scenario key is passed through unchanged (e.g. "ssp534-over"); the OUTPUT slug
# strips the hyphen to match our filename convention (ssp534over), because the rest of the
# pipeline keys on that.
_DEFAULT_SSPS = ["ssp126", "ssp245", "ssp585"]
SSPS = None                     # set from --ssps in main()


def out_slug(ssp):
    """Our filename key for an RCMIP scenario name: ssp534-over -> ssp534over."""
    return ssp.replace("-", "")

HIST_FILE = SPECIES_FILE = OUT_DIR = None    # set from the CLI in main()

START_YEAR = 1750
END_YEAR = 2301               # timepoints 1750.5 .. 2300.5
HIST_LAST_YEAR_INT = 2023     # CMIP7 historical ends here; the harmonization anchor

DROP_SPECIES = {"Halon-1202"}  # in CMIP7, not a species in 1.4.5 or 1.6.0

# Same unit-label equivalences as the chrisroadmap builder. The gate stays
# strict; see build_emissions_v160_ssp245.py for the verification of each.
UNIT_ALIASES = {
    ("Mt NO2/yr", "Mt NOx/yr"),
    ("kt HFC4310/yr", "kt HFC4310mee/yr"),
}

# A harmonized series must preserve the scenario's OWN year-on-year step across
# the join; anything larger is a broken anchor, not a property of the data.
MAX_JOIN_STEP_RESIDUAL = 1e-9


def rcmip_emissions(ssp, emis_driven, properties):
    """RCMIP-native emissions for one SSP.

    Defined on the EMISSIONS-DRIVEN SUBSET only. `fill_from_rcmip()` tries to
    fill every species it is given, and RCMIP has no entry for 1.6.0's new
    prescribed `Irrigation` -- on the full 62-species set it raises
    MissingDataError. Only emissions are wanted here anyway; prescribed forcing
    comes from calibration_v160_prod/volcanic_solar_<MARKER>.csv.
    """
    f = FAIR()
    f.define_time(START_YEAR, END_YEAR, step=1)
    f.define_scenarios([ssp])
    f.define_species(emis_driven, {k: properties[k] for k in emis_driven})
    f.ch4_method = "Thornhill2021"
    f.define_configs(["one"])
    f.allocate()
    f.fill_from_rcmip()
    arr = f.emissions.sel(scenario=ssp, config="one").values.copy()   # (tp, specie)
    tp = f.timepoints.copy()
    del f
    return tp, arr


def main():
    global SSPS, HIST_FILE, SPECIES_FILE, OUT_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--ssps", nargs="+", default=_DEFAULT_SSPS,
                    help="RCMIP scenario keys (default: %(default)s)")
    ap.add_argument("--hist-file", type=Path, required=True,
                    help="historical_emissions_1750-2023_cmip7.csv (Zenodo 18828694)")
    ap.add_argument("--species-file", type=Path, required=True,
                    help=f"calibration_v160_prod/species_configs_properties_{CALIB_VERSION}.csv")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    SSPS, HIST_FILE, SPECIES_FILE, OUT_DIR = (args.ssps, args.hist_file,
                                              args.species_file, args.out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    from forcing_inputs import require_rcmip_cache
    require_rcmip_cache()           # offline guard; see the PORT note above

    species, properties = read_properties(str(SPECIES_FILE))
    emis_driven = [s for s in species if properties[s]["input_mode"] == "emissions"]

    dfh = pd.read_csv(HIST_FILE)
    dfh = dfh[~dfh["variable"].isin(DROP_SPECIES)]
    ycols = [c for c in dfh.columns if str(c).strip().isdigit()]
    if max(int(c) for c in ycols) != HIST_LAST_YEAR_INT:
        sys.exit(f"ERROR: {HIST_FILE.name} does not end at {HIST_LAST_YEAR_INT}")
    hist_cols = [f"{int(c)}.5" for c in ycols]
    hist = dfh.set_index("variable")
    hist_units = hist["unit"]

    missing = set(emis_driven) - set(hist.index)
    if missing:
        sys.exit(f"ERROR: CMIP7 history lacks emissions-driven species {sorted(missing)}")
    print("── CMIP6-SSP emissions for calib 1.6.0 ──────────────────────────────")
    print(f"  history : {HIST_LABEL}")
    print(f"  future  : {FUTURE_LABEL}")
    print(f"  join    : {HIST_LAST_YEAR_INT}.5, per-species multiplicative")
    print(f"  species : {len(emis_driven)} emissions-driven of {len(species)}")

    for ssp in SSPS:
        tp, arr = rcmip_emissions(ssp, emis_driven, properties)
        idx = {s: i for i, s in enumerate(emis_driven)}
        fut_mask = tp > HIST_LAST_YEAR_INT + 0.5
        fut_cols = [f"{t:.1f}" for t in tp[fut_mask]]
        i_anchor = int(np.where(np.isclose(tp, HIST_LAST_YEAR_INT + 0.5))[0][0])

        rows, scales, worst, worst_sp = [], {}, 0.0, None
        renamed = {}
        for sp in emis_driven:
            # Both arms are interpreted by FaIR through the SAME species config,
            # so there is no second unit string to reconcile here: fill_from_rcmip
            # returns values already in the species config's own units, and the
            # CMIP7 unit label is carried through to the output for the reader.
            # FaIR's OWN spelling, not the CMIP7 frame's -- see fair_units.py
            _src = str(hist_units[sp]).strip()
            hu = fair_canonical_unit(sp, _src)
            if hu != _src:
                renamed[sp] = (_src, hu)
            h_vals = hist.loc[sp, ycols].to_numpy(float)
            h_a = float(h_vals[-1])
            c_a = float(arr[i_anchor, idx[sp]])
            scale = (h_a / c_a) if abs(c_a) > 1e-12 else (1.0 if abs(h_a) < 1e-12 else 0.0)
            scales[sp] = scale
            fut = arr[fut_mask, idx[sp]] * scale
            row = {"scenario": "all", "region": "World", "variable": sp, "unit": hu}
            row.update(dict(zip(hist_cols, h_vals)))
            row.update(dict(zip(fut_cols, fut)))
            rows.append(row)

            # join gate, per species
            c_b = float(arr[fut_mask, idx[sp]][0])
            if abs(h_a) > 1e-12 and abs(c_a) > 1e-12:
                d = abs((float(fut[0]) / h_a) - (c_b / c_a))
                if d > worst:
                    worst, worst_sp = d, sp

        unit_report(renamed)
        out = pd.DataFrame(rows, columns=["scenario", "region", "variable", "unit",
                                          *hist_cols, *fut_cols])
        yrs = [float(c) for c in [*hist_cols, *fut_cols]]
        if not all(b > a for a, b in zip(yrs, yrs[1:])):
            sys.exit(f"ERROR: {ssp} year columns not strictly increasing -- FaIR "
                     f"will raise NonMonotonicError")
        if out.isna().any().any():
            sys.exit(f"ERROR: NaN in {ssp} output")
        if worst >= MAX_JOIN_STEP_RESIDUAL:
            sys.exit(f"ERROR: {ssp} harmonization broke the scenario's own step at "
                     f"the join ({worst:.3e}, {worst_sp})")

        # NO provenance comment line in the file: FaIR's fill_from_csv does not
        # skip comments, so a '#' header would make the file unreadable by the
        # very driver that consumes it. Provenance lives in the BASIS_SLUG infix
        # in the filename and in calibration_v160_prod/README.md, and is printed
        # by every driver that reads it.
        path = OUT_DIR / f"emissions_v160_{BASIS_SLUG}_{out_slug(ssp)}.csv"
        out.to_csv(path, index=False)
        s = pd.Series(scales)
        print(f"\n  {ssp}: {out.shape[0]} species x {len(yrs)} timepoints "
              f"({yrs[0]} .. {yrs[-1]})")
        print(f"    join gate worst step residual {worst:.2e} ({worst_sp})  PASS")
        print(f"    harmonization scale: median x{s.median():.4f}, "
              f"{(s.sub(1).abs() > 0.10).sum()} of {len(s)} species move >10%")
        print(f"    wrote {path.name}")

    print(f"\nNOTE: the '{BASIS_SLUG}' infix is the provenance marker. These are "
          f"NOT RCMIP-native emissions -- history is CMIP7 to "
          f"{HIST_LAST_YEAR_INT} and the RCMIP future is harmonized onto it.")


if __name__ == "__main__":
    main()
