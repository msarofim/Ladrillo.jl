"""
build_fair_mean_v160.py

Ensemble-mean GMST and OHC trajectories for the SLR-RFF-BRICK (Ladrillo) side,
on the PRODUCTION fair-calibrate **1.6.0** stack. The 1.6.0 analogue of
`build_fair_mean_v145.py`, deliberately identical in output schema, baselining
and aggregation so the two are interchangeable for the consumer:

  fair_mean_gmst[_TAG].csv   year, gmst_C      (anomaly vs 1850-1900, K)
  fair_mean_ohc[_TAG].csv    year, ohc_1e22J   (anomaly vs 1850-1900, 1e22 J)

Both are the MEAN over the 841 configs of each config's own 1850-1900-rebased
trajectory -- rebase first, then average, exactly as the v145 driver does.
`stochastic_run` is left at whatever the posterior sets (both 1.4.5 and 1.6.0
set it True with per-config seeds), again matching v145, so the before/after
comparison the Ladrillo session will make is like-for-like.

WHAT DIFFERS FROM THE v145 DRIVER, AND WHY
  1. `--calib-dir` defaults to `calibration_v160_prod/`, the production stack:
     stock 1.6.0 Land use (`input_mode=forcing`), not the `calculated` deviation
     that `calibration_v160/` carries for the van-Vuuren-vs-RFF comparison.
  2. Prescribed forcing comes from a PER-MARKER file, `volcanic_solar_<MARKER>.csv`,
     carrying Volcanic + Solar + Land use + Irrigation. The marker must be named
     explicitly; see `--marker` below.
  3. The prescribed-forcing scaling is a LOOP over every `input_mode=forcing`
     species, and reports each one. The v145 driver hardcodes Volcanic and Solar,
     which silently leaves 1.6.0's new Irrigation UNSCALED -- FaIR applies
     `forcing_scale` internally only to CALCULATED species. A missing scale looks
     exactly like a genuinely cooler posterior, which is the failure mode that
     cost this project a week; commit fbb1375 made the same shape loud in the
     ERF drivers.

THE MARKER. 1.6.0 ships land-use and irrigation forcing for the seven CMIP7
ScenarioMIP markers; our drivers run CMIP6 SSPs and no CMIP7 marker IS a CMIP6
SSP. The two series are bit-identical across all seven markers through 2023 and
diverge only after (0.0895 and 0.0342 W/m2 by 2100), so the choice cannot touch
a historical result. `ssp245harm` uses **M**: chrisroadmap's file is
MESSAGE-GLOBIOM SSP2-4.5 and the CMIP7 M marker is the SSP2-based medium-forcing
marker, a direct correspondence rather than a stretch. ssp126/ssp585 have no
such clean counterpart and are an open question.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_fair_mean_v160.py @ 8984a13. The numerical code is copied
UNCHANGED; only the I/O plumbing differs:
  * --calib-dir and --out-dir are REQUIRED and resolve against the current
    directory (the original resolved them against the FaIRtoFrEDI root and
    defaulted --out-dir to SLR-RFF-BRICK/data/observations).
  * --emissions-file resolves against the current directory.
Normally driven by build_forcing.py. By hand (the ssp245harm calibration forcing
-- note the MIX of two directories, which is the documented trap):
    python build_fair_mean_v160.py \\
        --calib-dir WORK/calibration_v160_prod \\
        --emissions-file WORK/calibration_v160/emissions_v160_cmip7_ssp245_harmonized.csv \\
        --tag ssp245harm --scenario-label ssp245_harmonized --marker M \\
        --out-dir OUT
"""

from pathlib import Path
import argparse
import sys

import numpy as np
import pandas as pd
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fair_basis import declare, require_basis

# ── Named constants: every label, filename and printed line derives from these ─
START_YEAR = 1750
END_YEAR = 2301
BASELINE_START = 1850
BASELINE_END = 1900

MODEL_LABEL = "FaIR 2.2.4"
CALIB_LABEL = "calib 1.6.0"
# ⚠ HIST_LABEL WAS A CONSTANT HERE, AND THAT WAS THE BUG (fixed 2026-10-02).
# This driver takes --emissions-file, so the history basis is whatever the caller
# hands it -- but the provenance string asserted "CMIP7 historical 1750-2023"
# unconditionally, so pointing it at a v145 or a van Vuuren Smith-history file
# produced an output labelled CMIP7 that was not. The label is now DERIVED from
# the file by fair_basis.classify_basis, and the incoherent pair is refused.

# Sanity band for 2024 GMST anomaly vs 1850-1900. IGCC 2024 puts observed
# warming at ~1.24 K over 2015-2024 and ~1.5 K in 2024 itself; a driver that
# lands outside this has a broken forcing or emissions input, not a new result.
GMST_2024_LO, GMST_2024_HI = 0.9, 2.0


def resolve_calib(calib_dir, marker):
    d = Path(calib_dir).resolve()

    def one(pat):
        hits = sorted(d.glob(pat))
        if len(hits) != 1:
            sys.exit(f"ERROR: expected exactly one {pat} in {d}, found {len(hits)}")
        return hits[0]

    params = one("calibrated_constrained_parameters_*.csv")
    species = one("species_configs_properties_*.csv")
    # ⚠ MARKER-FREE MODE (2026-09-02). `--marker none` selects an unsuffixed
    # volcanic_solar.csv, which is how calibration_v160/ ships: Land use is
    # `calculated` from the scenario's own cumulative CO2 AFOLU, and Irrigation is
    # ONE shared trajectory (the 7-marker mean). That removes the CMIP7 marker
    # choice from the SSP runs entirely.
    # ⚠⚠ IT USES THE STOCK POSTERIOR OFF-DESIGN. The constrained parameters are
    # BYTE-IDENTICAL to the prescribed-land-use calibration (verified), so 1.6.0's
    # posterior was never constrained under `calculated` land use. That is only
    # safe if `calculated` reproduces the prescribed series over the CONSTRAINING
    # period; check it before trusting an absolute number, and see
    # calibration_v160/README.md, which states the deviation and its scope.
    forcing = (d / "volcanic_solar.csv") if str(marker).lower() in ("none", "shared", "") \
        else d / f"volcanic_solar_{marker}.csv"
    if not forcing.exists():
        avail = sorted(p.stem.replace("volcanic_solar_", "")
                       for p in d.glob("volcanic_solar_*.csv"))
        sys.exit(f"ERROR: no {forcing.name} in {d}. Available markers: {avail}")
    version = params.stem.replace("calibrated_constrained_parameters_", "")
    return d, params, species, forcing, version


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calib-dir", required=True,
                    help="calibration directory (the production stack is "
                         "calibration_v160_prod, built by "
                         "setup_calibration_v160_production.py)")
    ap.add_argument("--emissions-file", type=Path, required=True)
    ap.add_argument("--marker", required=True,
                    help="CMIP7 ScenarioMIP marker supplying prescribed land use "
                         "and irrigation (VL LN L ML M H HL). ssp245harm uses M.")
    ap.add_argument("--tag", default="",
                    help="output filename infix: fair_mean_gmst[_TAG].csv")
    ap.add_argument("--scenario-label", required=True,
                    help="FaIR scenario dimension label")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--allow-legacy-basis", action="store_true",
                    help="permit a non-CMIP7 emissions file under calib 1.6.0 -- "
                         "the measured input mismatch. The output's provenance "
                         "string is marked WAIVED when this is used.")
    args = ap.parse_args()

    CAL_DIR, PARAMS_FILE, SPECIES_FILE, FORCING_FILE, CAL_VER = resolve_calib(
        args.calib_dir, args.marker)
    emis = args.emissions_file
    if not emis.exists():
        sys.exit(f"ERROR: missing emissions file {emis}")
    suffix = f"_{args.tag}" if args.tag else ""
    out_gmst = args.out_dir / f"fair_mean_gmst{suffix}.csv"
    out_ohc = args.out_dir / f"fair_mean_ohc{suffix}.csv"

    basis = require_basis(emis, calib_version=CAL_VER,
                          allow_legacy=args.allow_legacy_basis)
    provenance = (f"{MODEL_LABEL} ({CALIB_LABEL}); {basis.text}; "
                  f"marker {args.marker}; "
                  f"anomaly vs {BASELINE_START}-{BASELINE_END}; "
                  f"mean over configs")

    print(f"build_fair_mean_v160.py — {MODEL_LABEL} ({CALIB_LABEL}) ensemble means")
    print(f"  calibration : {CAL_VER}  ({CAL_DIR.name})")
    print(f"  emissions   : {emis.name}")
    print(f"  basis       : {basis.text}")
    print(f"  forcing     : {FORCING_FILE.name}  (marker {args.marker})")
    print(f"  scenario    : {args.scenario_label}  -> tag '{args.tag or '(none)'}'")

    df_configs = pd.read_csv(PARAMS_FILE, index_col=0)
    configs = df_configs.index.tolist()
    species, properties = read_properties(str(SPECIES_FILE))
    prescribed = [s for s in species if properties[s]["input_mode"] == "forcing"]
    print(f"  {len(configs)} configs, {len(species)} species, "
          f"{len(prescribed)} prescribed-forcing: {prescribed}")
    # The prescribed/emissions split in one printable line, and gated: a
    # prescribed specie with no forcing_scale column would run UNSCALED, which
    # looks exactly like a cooler posterior. The prescribed SET changes with the
    # calibration DIRECTORY, not just the version (calibration_v160 runs Land use
    # `calculated`, calibration_v160_prod prescribes it per marker), so the
    # declaration goes into the provenance string, not only into stdout.
    decl = declare(species, properties, FORCING_FILE, df_configs)
    print(f"  {decl.text}")
    provenance = f"{provenance}; {decl.text}"

    f = FAIR()
    f.define_time(START_YEAR, END_YEAR, step=1)
    f.define_scenarios([args.scenario_label])
    f.define_species(species, properties)
    f.ch4_method = "Thornhill2021"
    f.define_configs(configs)
    f.allocate()

    # ⚠ EACH CSV IS READ ON ITS OWN SCENARIO KEY, which until 2026-10-02 was
    # hardcoded to "all" for both. `fill_from_csv` matches the `scenario` column
    # and fills NOTHING when it does not match, then dies in the unit lookup with
    # `IndexError: index 0 is out of bounds` -- an error that names a unit and
    # says nothing about scenarios. The v160 SSP and RFF-SP payloads do carry
    # "all", but the van Vuuren files carry their own prose name ("Medium - SSP2
    # (Marker)"), which the vv drivers read to identify the marker. So the key is
    # now read from the file instead of assumed.
    def _scen_key(path, col):
        hdr = pd.read_csv(path, nrows=1)
        hit = [c for c in hdr.columns if c.strip().lower() == col]
        if not hit:
            sys.exit(f"ERROR: {Path(path).name} has no {col!r} column")
        return str(hdr[hit[0]].iloc[0])

    scen_e = _scen_key(emis, "scenario")
    scen_f = _scen_key(FORCING_FILE, "scenario")

    def _load(path, kind, key):
        t = FAIR()
        t.define_time(START_YEAR, END_YEAR, step=1)
        t.define_scenarios([key])
        t.define_species(species, properties)
        t.ch4_method = "Thornhill2021"
        t.define_configs(configs)
        t.allocate()
        t.fill_from_csv(**{f"{kind}_file": str(path),
                           ("forcing" if kind == "emissions" else "emissions") + "_file": None})
        arr = getattr(t, kind).sel(scenario=key).values.copy()
        del t
        return arr

    emis_arr = _load(emis, "emissions", scen_e)
    forc_arr = _load(FORCING_FILE, "forcing", scen_f)
    print(f"  scenario keys read from the files: emissions {scen_e!r}, "
          f"forcing {scen_f!r}")

    # A prescribed specie the forcing CSV does not cover reads back all-NaN and
    # would poison the forcing sum. Fail rather than zero-fill: production must
    # not silently drop a forcing term.
    idx = {s: i for i, s in enumerate(species)}
    blank = [s for s in prescribed if np.isnan(forc_arr[:, :, idx[s]]).all()]
    if blank:
        sys.exit(f"ERROR: {FORCING_FILE.name} has no rows for prescribed species "
                 f"{blank}; production cannot run with a missing forcing term.")
    # Same for emissions-driven species.
    ed = [s for s in species if properties[s]["input_mode"] == "emissions"]
    miss = sorted({s for s in ed if np.isnan(emis_arr[:, :, idx[s]]).any()})
    if miss:
        sys.exit(f"ERROR: {emis.name} has NaN for emissions-driven species {miss}")

    f.emissions.loc[dict(scenario=args.scenario_label)] = emis_arr
    f.forcing.loc[dict(scenario=args.scenario_label)] = forc_arr

    f.fill_species_configs(str(SPECIES_FILE))
    f.override_defaults(str(PARAMS_FILE))

    # FaIR applies forcing_scale internally to CALCULATED species only. Every
    # prescribed one must be scaled here, and reported either way -- a silent
    # skip is indistinguishable from a cooler posterior.
    for sp in prescribed:
        col = f"forcing_scale[{sp}]"
        if col in df_configs.columns:
            fill(f.forcing, f.forcing.sel(specie=sp)
                 * df_configs[col].values.squeeze(), specie=sp)
            print(f"  scaled prescribed forcing: {sp}")
        else:
            sys.exit(f"ERROR: prescribed specie {sp!r} has no {col} in "
                     f"{PARAMS_FILE.name}; it would run UNSCALED.")

    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    initialise(f.ocean_heat_content_change, 0)

    print(f"\nRunning FaIR ({len(configs)} configs, single scenario) ...")
    f.run(progress=False)
    print("  done.")

    years = f.timebounds.astype(int)
    pi = (years >= BASELINE_START) & (years <= BASELINE_END)

    temp_raw = f.temperature.sel(layer=0).values[:, 0, :]      # (T, cfg)
    gmst_mean = (temp_raw - temp_raw[pi, :].mean(axis=0)[None, :]).mean(axis=1)

    ohc_raw = f.ocean_heat_content_change.values[:, 0, :]      # (T, cfg), J
    ohc_mean = ((ohc_raw - ohc_raw[pi, :].mean(axis=0)[None, :]) / 1e22).mean(axis=1)

    i2024 = int(np.where(years == 2024)[0][0])
    g2024 = float(gmst_mean[i2024])
    ok = GMST_2024_LO < g2024 < GMST_2024_HI
    print(f"\n  gate: GMST 2024 = {g2024:+.4f} K vs "
          f"{BASELINE_START}-{BASELINE_END}  "
          f"(expect {GMST_2024_LO}..{GMST_2024_HI})  {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("ERROR: 2024 GMST outside the plausible band -- suspect the "
                 "emissions or prescribed-forcing input, not a new result.")

    print("\nGMST mean (K, rel 1850-1900) / OHC mean (1e22 J):")
    for y in [1850, 1900, 1971, 2000, 2024, 2050, 2100, 2300]:
        i = int(np.where(years == y)[0][0])
        print(f"  {y}: {gmst_mean[i]:+.4f} K   {ohc_mean[i]:+9.3f}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"year": years, "gmst_C": gmst_mean}).to_csv(out_gmst, index=False)
    pd.DataFrame({"year": years, "ohc_1e22J": ohc_mean}).to_csv(out_ohc, index=False)
    print(f"\nWrote {out_gmst}")
    print(f"Wrote {out_ohc}")
    print(f"\nPROVENANCE: {provenance}")


if __name__ == "__main__":
    main()
