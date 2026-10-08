"""
build_fair_cube_vv_v160.py

PER-CONFIG GMST/OHC cubes and their ensemble means for one van Vuuren
ScenarioMIP-CMIP7 MARKER, on the PRODUCTION fair-calibrate **1.6.0** stack.

The van Vuuren counterpart of `build_fair_cube_v160.py`. Writes the SAME schema
into SLR-RFF-BRICK/data/observations/, under a `vv<MARKER>` scenario key:
    fair_cube_gmst_vv<MARKER>_raw.csv    year, cfg_<id>...   K rel 1850-1900
    fair_cube_ohc_vv<MARKER>_raw.csv     year, cfg_<id>...   1e22 J rel 1850-1900
    fair_mean_gmst_vv<MARKER>.csv        year, gmst_C
    fair_mean_ohc_vv<MARKER>.csv         year, ohc_1e22J

Schema, year range (1850..2300 for cubes, full range for means), 6-decimal cube
rounding, 1850-1900 baselining and rebase-then-average order all match the SSP
cubes exactly, so `julia/scope_slr_fair_uncertainty.jl --ssp=vvH` and every other
consumer needs no edit -- they already resolve their inputs by scenario key.

WHY THE MARKER FORCING IS THE POINT
fair-calibrate 1.6.0 ships prescribed land-use and irrigation forcing for the
SEVEN CMIP7 ScenarioMIP markers. `calibration_v160_prod/README.md` states the
problem this creates for the SSP cubes plainly: "no CMIP7 marker **is** a CMIP6
SSP, so a driver must name the marker it uses ... It does affect projections to
2300 and is an open question for ssp126/ssp585." The SSP cubes therefore run on
an APPROXIMATE mapping (ssp126->L, ssp245->M, ssp585->H) whose ambiguity was
measured at 0.022-0.094 K at 2100/2300. Here that mapping is EXACT: each van
Vuuren marker uses its own `volcanic_solar_<MARKER>.csv`, and the ambiguity is
identically zero. This is the one thing the van Vuuren set gets for free that
the SSP set cannot.

⚠ NOT `calibration_v160/`. That directory deliberately deviates from stock 1.6.0
(Land use forced to `calculated`; Irrigation as one shared 7-marker mean) for the
van-Vuuren-vs-RFF PERCENTILE comparison, and is wrong for production -- the
land-use deviation alone is 0.025 K of historical GMST fit. The existing
`fair_outputs/vanvuuren_ext_harmonized_erf_cube_v160.npz` was built against it
and against the single generic `volcanic_solar.csv`, so it is NOT the ERF
counterpart of these cubes.

STOCHASTIC IS OFF BY DEFAULT, AND THAT IS A DEPARTURE FROM THE SSP CUBES.
The 1.6.0 posterior sets `stochastic_run=True` with 841 distinct seeds, so the
SSP cubes carry FaIR's internal-variability noise. For a cube whose only job is
to drive an SLR model, that noise is not an uncertainty about the forced
response: it inflates the joint band with a quantity that cannot be separated
from parameter spread afterwards, and thermal expansion INTEGRATES it, so it
accumulates rather than averaging out over 450 years. `--stochastic=on`
reproduces the SSP convention; `scripts/diag_vv_stochastic_ab.py` is the
measurement that decides which is right, and its numbers belong in the
CHANGELOG entry for whichever is chosen.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_fair_cube_vv_v160.py @ 8984a13. The numerical code is copied
UNCHANGED (including --stochastic defaulting to OFF); only the I/O plumbing
differs:
  * --calib-dir and --out-dir are REQUIRED and resolve against the current
    directory (the original hard-coded FaIRtoFrEDI/calibration_v160_prod and
    SLR-RFF-BRICK/data/observations).
  * the CMIP7-basis emissions are read from --emissions-dir (default: the
    calibration dir, as in the original); the legacy Smith-basis directory is
    --legacy-emissions-dir and is only read with --legacy-basis.
Normally driven by build_forcing.py. By hand:
    python build_fair_cube_vv_v160.py --marker H \
        --calib-dir WORK/calibration_v160_prod --out-dir OUT
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

# ── Named constants ──────────────────────────────────────────────────────────
START_YEAR = 1750
END_YEAR = 2301
BASELINE_START, BASELINE_END = 1850, 1900
CUBE_FIRST_YEAR = 1850
CUBE_LAST_YEAR = 2300
CUBE_DECIMALS = 6

CALIB_VERSION = "1.6.0"
MODEL_LABEL = "FaIR 2.2.4"
SCEN_PREFIX = "vv"              # scenario key = SCEN_PREFIX + MARKER, e.g. "vvH"

CAL_DIR = None                  # set from --calib-dir in main()

# ── THE EMISSIONS BASIS IS NOW A CHOICE, AND THE DEFAULT CHANGED 2026-10-02 ──
# This driver used to read data/vanvuuren/spliced_ext_harmonized/ unconditionally
# while its provenance string asserted "CMIP7 historical 1750-2023". That claim
# was FALSE: those files were built for calib 1.4.5 and their 1750.5-2023.5 block
# is SMITH 2024 history, so the 1.6.0 posterior was running on v1.4.5 history --
# the rung-R1 input mismatch. Measured on markers M and H, holding everything
# else: GMST 2015-2024 bias vs IGCC -0.1092 -> -0.0608 K and RMSE 0.1317 -> 0.1217 K.
# ⚠ The mechanism is a 0.065 K COOLER 1850-1900 reference, not extra warming:
# absolute 2100 GMST moves -0.010 K and absolute total ERF <= 0.022 W/m2, so no
# absolute-threshold statistic changes (P(ERF>8.5) @2150 marker H is 42.212 % on
# both bases). An earlier version of this comment claimed +0.116 W/m2 and a
# 6.4 pp exceedance move; both are RETRACTED -- they re-referenced ERF to
# 1850-1900 before comparing it to an ABSOLUTE 8.5 W/m2 threshold.
# The default is now the CMIP7 basis; the old basis is still reachable, but only
# by asking for it.
# ⚠ Any vv cube built before 2026-10-02 is on the LEGACY basis.
CMIP7_EMIS_GLOB = "emissions_v160_cmip7harm_vv*.csv"   # dir: --emissions-dir
LEGACY_EMIS_GLOB = "*.csv"      # dir: --legacy-emissions-dir (spliced_ext_harmonized)

# The van Vuuren scenario NAME carried in each spliced CSV's own `scenario`
# column -> the CMIP7 marker code that names its forcing file. Declared here
# rather than parsed out of the prose name: "High-to-Low - SSP5 (Marker)" -> HL
# is a reading, and a reading belongs in a constant a human can check, not in a
# string operation that will one day quietly map a new scenario to the wrong
# forcing. Gated below against BOTH the emissions dir and the calibration dir.
SCENARIO_TO_MARKER = {
    "High - SSP3 (Marker)":            "H",
    "High-to-Low - SSP5 (Marker)":     "HL",
    "Medium - SSP2 (Marker)":          "M",
    "Medium-to-Low - SSP2 (Marker)":   "ML",
    "Low - SSP2 (Marker)":             "L",
    "Low-to-Negative - SSP2 (Marker)": "LN",
    "Very Low - SSP1 (Marker)":        "VL",
}

MEAN_MATCH_TOL = 1e-12          # see build_fair_cube_v160.py: float64 repr floor, not picked

# First-principles band on 2024 GMST vs 1850-1900 (IGCC 2024 ~1.24 K over
# 2015-2024). Outside this, suspect the input, not the result.
GMST_2024_LO, GMST_2024_HI = 0.9, 2.0


def resolve_marker(marker, emis_dir, emis_glob):
    """Map a marker code back to its emissions file, gating the whole table.

    Every check here is against something on disk. The mapping is a BIJECTION
    over exactly the seven markers, each marker has a forcing file, and each
    scenario name is actually present in the emissions directory -- so a typo in
    SCENARIO_TO_MARKER is an error at startup, not a cube built on the wrong
    scenario's land use. `emis_glob` is needed because the CMIP7-basis files live
    in the calibration directory alongside the SSP payloads, whose `scenario`
    column is "all" and is not a marker.
    """
    markers = sorted(SCENARIO_TO_MARKER.values())
    if len(set(markers)) != len(SCENARIO_TO_MARKER):
        sys.exit(f"ERROR: SCENARIO_TO_MARKER is not a bijection: {markers}")

    csvs = sorted(emis_dir.glob(emis_glob))
    if not csvs:
        sys.exit(f"ERROR: no emissions CSVs matching {emis_glob} in {emis_dir}. "
                 f"The CMIP7-basis van Vuuren files are built by "
                 f"scripts/build_emissions_v160_cmip7harm_vv.py.")
    scen_by_path = {p: pd.read_csv(p, nrows=1)["scenario"].iloc[0] for p in csvs}

    unknown = sorted(set(scen_by_path.values()) - set(SCENARIO_TO_MARKER))
    if unknown:
        sys.exit(f"ERROR: {emis_dir.name}/{emis_glob} holds scenario(s) absent "
                 f"from SCENARIO_TO_MARKER: {unknown}")
    missing = sorted(set(SCENARIO_TO_MARKER) - set(scen_by_path.values()))
    if missing:
        sys.exit(f"ERROR: SCENARIO_TO_MARKER names scenario(s) with no emissions "
                 f"file in {emis_dir.name}/{emis_glob}: {missing}")

    for m in markers:
        fp = CAL_DIR / f"volcanic_solar_{m}.csv"
        if not fp.exists():
            sys.exit(f"ERROR: marker {m} has no {fp.name} in {CAL_DIR.name}")

    hit = [(p, s) for p, s in scen_by_path.items()
           if SCENARIO_TO_MARKER[s] == marker]
    if len(hit) != 1:
        sys.exit(f"ERROR: marker {marker!r} matched {len(hit)} emissions files")
    return hit[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--marker", required=True,
                    choices=sorted(set(SCENARIO_TO_MARKER.values())))
    ap.add_argument("--stochastic", default="off", choices=["on", "off"],
                    help="FaIR internal variability. Default OFF -- see module "
                         "docstring; 'on' reproduces the SSP-cube convention.")
    ap.add_argument("--calib-dir", type=Path, required=True,
                    help="the production stack, calibration_v160_prod")
    ap.add_argument("--emissions-dir", type=Path, default=None,
                    help="where emissions_v160_cmip7harm_vv<MARKER>.csv live "
                         "(default: --calib-dir)")
    ap.add_argument("--legacy-emissions-dir", type=Path, default=None,
                    help="spliced_ext_harmonized/ (SMITH 2024 history); read "
                         "only with --legacy-basis")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--legacy-basis", action="store_true",
                    help="read data/vanvuuren/spliced_ext_harmonized/ (SMITH 2024 "
                         "history) instead of the CMIP7-basis files. This pairs "
                         "the 1.6.0 posterior with v1.4.5 history -- the measured "
                         "input mismatch. Use ONLY to reproduce a pre-2026-10-02 "
                         "cube, and say so wherever the number is used.")
    args = ap.parse_args()
    global CAL_DIR
    CAL_DIR = args.calib_dir
    CMIP7_EMIS_DIR = args.emissions_dir or CAL_DIR
    LEGACY_EMIS_DIR = args.legacy_emissions_dir
    if args.legacy_basis and LEGACY_EMIS_DIR is None:
        sys.exit("ERROR: --legacy-basis needs --legacy-emissions-dir")

    marker = args.marker
    scen_key = f"{SCEN_PREFIX}{marker}"
    emis_dir, emis_glob = ((LEGACY_EMIS_DIR, LEGACY_EMIS_GLOB) if args.legacy_basis
                           else (CMIP7_EMIS_DIR, CMIP7_EMIS_GLOB))
    emis, scen_name = resolve_marker(marker, emis_dir, emis_glob)
    # Refuses the (CMIP7-conditioned calibration, non-CMIP7 emissions) pair
    # unless --legacy-basis was passed, and marks the Basis WAIVED when it was,
    # so the waiver reaches the provenance string instead of dying at the flag.
    basis = require_basis(emis, calib_version=CALIB_VERSION,
                          allow_legacy=args.legacy_basis)

    params = CAL_DIR / f"calibrated_constrained_parameters_{CALIB_VERSION}.csv"
    species_file = CAL_DIR / f"species_configs_properties_{CALIB_VERSION}.csv"
    forcing_file = CAL_DIR / f"volcanic_solar_{marker}.csv"
    for p in (params, species_file, forcing_file, emis):
        if not p.exists():
            sys.exit(f"ERROR: missing {p}")

    stoch = args.stochastic == "on"
    df_configs = pd.read_csv(params, index_col=0)
    configs = df_configs.index.tolist()
    species, properties = read_properties(str(species_file))
    prescribed = [s for s in species if properties[s]["input_mode"] == "forcing"]

    # ⚠ THE HISTORY LABEL IS DERIVED FROM THE FILE, NOT ASSERTED. It used to be
    # the literal "CMIP7 historical 1750-2023" while the driver read a Smith-2024
    # file -- a provenance string that was simply false. `basis.text` is computed
    # from the emissions filename by fair_basis.classify_basis, so changing the
    # input changes the label by construction.
    decl = declare(species, properties, forcing_file, df_configs)
    provenance = (
        f"{MODEL_LABEL} (calib {CALIB_VERSION}); van Vuuren ScenarioMIP-CMIP7 "
        f"{scen_name}; {basis.text}; published post-2100 extension "
        f"(Zenodo 20713982 v1.1.1), truncated at 2300; prescribed "
        f"land use + irrigation from its OWN CMIP7 marker {marker} (exact, not "
        f"mapped); {decl.text}; stochastic_run={stoch}; anomaly vs "
        f"{BASELINE_START}-{BASELINE_END}")

    print(f"── per-config cube, {scen_key} ──────────────────────────────────────")
    print(f"  calibration : {CALIB_VERSION}  ({CAL_DIR.name})")
    print(f"  scenario    : {scen_name}")
    print(f"  emissions   : {emis.name}")
    print(f"  basis       : {basis.text}")
    print(f"  forcing     : {forcing_file.name}  (marker {marker}, EXACT)")
    print(f"  {decl.text}")
    print(f"  stochastic  : {'ON' if stoch else 'OFF'}")
    print(f"  {len(configs)} configs, {len(species)} species, "
          f"prescribed: {prescribed}")

    f = FAIR()
    f.define_time(START_YEAR, END_YEAR, step=1)
    f.define_scenarios([scen_key])
    f.define_species(species, properties)
    f.ch4_method = "Thornhill2021"
    f.define_configs(configs)
    f.allocate()

    # volcanic_solar_<MARKER>.csv carries scenario="all"; the emissions CSV
    # carries the van Vuuren prose name. Neither is scen_key, so both are loaded
    # through a temp FaIR on their own scenario axis and copied in -- the same
    # load-and-broadcast the SSP builder and the VV ERF driver both use, because
    # passing a mismatched scenario name to fill_from_csv gives a 0-row lookup.
    t = FAIR()
    t.define_time(START_YEAR, END_YEAR, step=1)
    t.define_scenarios(["all"])
    t.define_species(species, properties)
    t.ch4_method = "Thornhill2021"
    t.define_configs(configs)
    t.allocate()
    t.fill_from_csv(emissions_file=None, forcing_file=str(forcing_file))
    forc_arr = t.forcing.sel(scenario="all").values
    del t

    e = FAIR()
    e.define_time(START_YEAR, END_YEAR, step=1)
    e.define_scenarios([scen_name])
    e.define_species(species, properties)
    e.ch4_method = "Thornhill2021"
    e.define_configs(configs)
    e.allocate()
    e.fill_from_csv(emissions_file=str(emis), forcing_file=None)
    emis_arr = e.emissions.sel(scenario=scen_name).values
    del e

    idx = {s: i for i, s in enumerate(species)}
    blank = [s for s in prescribed if np.isnan(forc_arr[:, :, idx[s]]).all()]
    if blank:
        sys.exit(f"ERROR: {forcing_file.name} has no rows for prescribed {blank}")
    ed = [s for s in species if properties[s]["input_mode"] == "emissions"]
    miss = sorted({s for s in ed if np.isnan(emis_arr[:, :, idx[s]]).any()})
    if miss:
        sys.exit(f"ERROR: {emis.name} has NaN for emissions-driven {miss}")

    f.emissions.loc[dict(scenario=scen_key)] = emis_arr
    f.forcing.loc[dict(scenario=scen_key)] = forc_arr
    f.fill_species_configs(str(species_file))
    f.override_defaults(str(params))

    # override_defaults sets stochastic_run from the posterior (True for 1.6.0).
    # Set it AFTER, and REPORT the value actually in force -- a run whose noise
    # setting is not printed cannot be told apart from one that silently kept
    # the posterior's.
    f.climate_configs["stochastic_run"][:] = stoch
    got = bool(np.all(f.climate_configs["stochastic_run"].values == stoch))
    if not got:
        sys.exit("ERROR: stochastic_run did not take the requested value")
    print(f"  climate_configs['stochastic_run'] = {stoch} (all {len(configs)} configs)")

    for sp in prescribed:
        col = f"forcing_scale[{sp}]"
        if col not in df_configs.columns:
            sys.exit(f"ERROR: prescribed {sp!r} has no {col}; would run UNSCALED")
        fill(f.forcing, f.forcing.sel(specie=sp)
             * df_configs[col].values.squeeze(), specie=sp)
        print(f"  scaled prescribed forcing: {sp}")

    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    initialise(f.ocean_heat_content_change, 0)

    print(f"\nRunning FaIR ({len(configs)} configs) ...")
    f.run(progress=False)

    years = f.timebounds.astype(int)
    pi = (years >= BASELINE_START) & (years <= BASELINE_END)
    gmst = f.temperature.sel(layer=0).values[:, 0, :]
    gmst = gmst - gmst[pi, :].mean(axis=0)[None, :]
    ohc = f.ocean_heat_content_change.values[:, 0, :]
    ohc = (ohc - ohc[pi, :].mean(axis=0)[None, :]) / 1e22

    i2024 = int(np.where(years == 2024)[0][0])
    g2024 = float(gmst[i2024].mean())
    ok = GMST_2024_LO < g2024 < GMST_2024_HI
    print(f"\n  gate: mean GMST 2024 = {g2024:+.4f} K   {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("ERROR: 2024 GMST outside the plausible band")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cy = (years >= CUBE_FIRST_YEAR) & (years <= CUBE_LAST_YEAR)
    written = {}
    for lab, arr, mean_col in [("gmst", gmst, "gmst_C"), ("ohc", ohc, "ohc_1e22J")]:
        cube = pd.DataFrame(np.round(arr[cy, :], CUBE_DECIMALS),
                            columns=[f"cfg_{c}" for c in configs])
        cube.insert(0, "year", years[cy])
        cp = args.out_dir / f"fair_cube_{lab}_{scen_key}_raw.csv"
        cube.to_csv(cp, index=False)
        mp = args.out_dir / f"fair_mean_{lab}_{scen_key}.csv"
        pd.DataFrame({"year": years, mean_col: arr.mean(axis=1)}).to_csv(mp, index=False)
        written[lab] = (cp, mp, arr)
        print(f"  wrote {cp.name}  ({cp.stat().st_size / 1e6:.1f} MB, "
              f"{len(configs)} cfg x {int(cy.sum())} yr)")
        print(f"  wrote {mp.name}")

    # [MEAN-MATCH] the mean file IS the cube's reduction, exactly, in memory.
    # ⚠ There is deliberately NO [CUBE-TRUNCATION] gate here: that one compares
    # against a quarantined pre-fix cube, and van Vuuren has no predecessor to
    # compare against. Saying so beats leaving a reader to wonder whether it was
    # dropped or forgotten.
    print("\n  [MEAN-MATCH] mean file vs cube reduction:")
    for lab, (cp, mp, arr) in written.items():
        col = "gmst_C" if lab == "gmst" else "ohc_1e22J"
        m_disk = pd.read_csv(mp).set_index("year")[col]
        d_mem = float(np.max(np.abs(m_disk.to_numpy() - arr.mean(axis=1))))
        cube_disk = pd.read_csv(cp).set_index("year")
        yy = cube_disk.index.to_numpy()
        d_disk = float(np.max(np.abs(cube_disk.to_numpy().mean(axis=1)
                                     - m_disk.loc[yy].to_numpy())))
        status = "PASS" if d_mem <= MEAN_MATCH_TOL else "FAIL"
        print(f"    {lab}: post-round-trip {d_mem:.3e} (<= {MEAN_MATCH_TOL:.0e}) "
              f"{status}   |   after {CUBE_DECIMALS}-dp cube rounding {d_disk:.3e}")
        if d_mem > MEAN_MATCH_TOL:
            sys.exit(f"ERROR: {lab} mean file is not the cube's own reduction")

    print(f"\nPROVENANCE: {provenance}")


if __name__ == "__main__":
    main()
