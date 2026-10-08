"""
build_fair_cube_v160.py

PER-CONFIG GMST/OHC cubes and their ensemble means for one CMIP6 SSP, on the
PRODUCTION fair-calibrate **1.6.0** stack. The 1.6.0 replacement for the cube
half of `run_fair_ssp_spread.py` plus `run_fair_ssps.py`.

Writes, into --out-dir (in the research repo: SLR-RFF-BRICK/data/observations/):
    fair_cube_gmst_<SSP>_raw.csv    year, cfg_<id>...   K rel 1850-1900
    fair_cube_ohc_<SSP>_raw.csv     year, cfg_<id>...   1e22 J rel 1850-1900
    fair_mean_gmst_<SSP>.csv        year, gmst_C
    fair_mean_ohc_<SSP>.csv         year, ohc_1e22J

Schema, year range (1850..2300 for cubes, full range for means), 6-decimal cube
rounding, 1850-1900 baselining and rebase-then-average order all match the v145
files exactly, so the Ladrillo consumers need no change.

WHAT IS DELIBERATELY NOT HERE
`run_fair_ssp_spread.py` also builds the **Coulon 2025 AIS arm** products --
`fair_coulon_arm_<ssp>.csv`, `fair_pctile_gmst_<ssp>.csv` and the hot-tail table.
Those are a controlled experiment whose configs are SELECTED by their 2300
warming against Coulon's +12.0/+14.5/+17.0 degC targets at the posterior-median
`ais_gmst_amp`. Re-running them under 1.6.0 means re-selecting the arm configs
against a moved ensemble and a moved L14 posterior -- an analysis decision, not a
mechanical port, and it depends on an L14 that does not exist yet.

CORRECTION (same day): an earlier version of this docstring also listed
`fair_cube_*_ssp585_spliced.csv` as a Coulon product. THAT WAS WRONG on two
counts. `run_fair_ssp_spread.py` does not write it at all -- it computes `spliced`
only to SELECT arms from (line 242 splices every config, line 246 selects) and
says so at its own cube block: "The SPLICE convention is applied in Julia ... not
baked in here." And the splice is separable from the selection: it is a pure
function of the raw cube and the ensemble-mean path, with no posterior and no
targets in it. `scripts/build_fair_cube_spliced_v160.py` emits it, gated by
reproducing the shipped v145 reference from v145 inputs first. Leaving it stale
would have failed a hard `@assert` in
`julia/scope_slr_fair_uncertainty.jl` ([SPLICE-MATCH], tol 2e-6) against a
0.38 K vintage gap -- correctly, but for the wrong reason.

EMISSIONS BASIS HAS CHANGED, AND NOT SILENTLY. The v145 cube driver took
emissions from `fill_from_rcmip()` -- RCMIP-native throughout, history included.
This driver reads `emissions_v160_cmip7harm_<SSP>.csv`: CMIP7 historical
1750-2023 with the RCMIP future harmonized onto it at 2023.5 (Marcus 2026-08-28).
The two are different objects. The `cmip7harm` infix and the provenance line
printed on every run are the marker.

MARKER MAPPING (Marcus 2026-08-28, by 2100 forcing level):
    ssp126 -> L    ssp245 -> M    ssp585 -> H
Measured cost of the ambiguity, emissions held fixed: 0.000 K through 2023, then
at 2100/2300 ssp126 0.022/0.045, ssp245 0.017/0.059, ssp585 0.037/0.094 K.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_fair_cube_v160.py @ 8984a13. The numerical code is copied
UNCHANGED; only the I/O plumbing differs:
  * --calib-dir and --out-dir are REQUIRED and resolve against the current
    directory (the original defaulted to FaIRtoFrEDI/calibration_v160_prod and
    SLR-RFF-BRICK/data/observations).
  * --emissions-dir resolves against the current directory.
  * the [CUBE-TRUNCATION] gate's reference directory (a research-repo
    quarantine) is now --truncation-reference-dir; without it the gate is
    skipped AND SAYS SO. It never touched the arithmetic.
Normally driven by build_forcing.py. By hand:
    python build_fair_cube_v160.py --ssp ssp245 \
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

# ── Named constants ──────────────────────────────────────────────────────────
START_YEAR = 1750
END_YEAR = 2301
BASELINE_START, BASELINE_END = 1850, 1900
CUBE_FIRST_YEAR = 1850          # v145 cubes start here, not at START_YEAR
# ...AND they STOP here, one year short of END_YEAR. run_fair_ssp_spread.py:103
# had BOTH bounds (`CUBE_Y0, CUBE_Y1 = 1850, 2300`) and sliced two-sided at :260.
# The first port carried only the lower one, so the v160 cubes shipped 452 rows
# where every downstream constant expects 451 -- caught by the hard 1850:2300
# assert at julia/scope_slr_fair_uncertainty.jl:186, which indexes the cube
# POSITIONALLY and would otherwise have shifted a projection band by a year.
# 2301 is a legitimate FaIR timebound; it was simply never part of this
# convention. The `fair_mean_*` files ARE full-range in both vintages.
CUBE_LAST_YEAR = 2300
CUBE_DECIMALS = 6               # v145 cube rounding; keep so file sizes match

CALIB_VERSION = "1.6.0"
BASIS_SLUG = "cmip7harm"
MODEL_LABEL = "FaIR 2.2.4"

# The research repo kept the pre-fix (452-row) cubes in
# SLR-RFF-BRICK/outputs/quarantine/20260828_cube_year_2301; pass that directory
# as --truncation-reference-dir to re-run the [CUBE-TRUNCATION] gate.

# Marcus 2026-08-28: map each CMIP6 SSP to the CMIP7 ScenarioMIP marker at the
# same 2100 forcing level. ssp245harm already shipped on M; this keeps ssp245
# consistent with it.
## ⚠ ssp534over -> HL added 2026-09-02 (Marcus). SSP5-3.4-OS is SLEIP's overshoot scenario
## (egusphere-2026-3874) and HL is the only CMIP7 overshoot marker, so it matches the
## scenario's SHAPE; the alternative, H, would match its SSP5 socioeconomic family instead.
## The choice is worth 0.022-0.094 K at 2100/2300 on the other SSPs and is the documented
## "open question" class. It also keeps this arm on the SAME land-use basis as our vvHL van
## Vuuren marker, so the two overshoot comparisons are mutually consistent.
SSP_TO_MARKER = {"ssp126": "L", "ssp245": "M", "ssp585": "H", "ssp534over": "HL"}

# [MEAN-MATCH] the mean file and the cube must be the same reduction of the same
# run. The mean is compared AFTER a CSV round-trip, so the floor is float64
# REPRESENTATION noise, not zero: on a GMST of order 1-5 K one ulp is ~4e-16 and
# the observed residual is 2e-16 to 9e-16, i.e. 1-4 ulp. The tolerance is derived
# from that magnitude rather than picked -- 1e-12 is ~1000x the repr floor and
# still ~1e-12 of the quantity it guards, so it passes bit-noise and would catch
# any real reduction mismatch. It is NOT the disk-rounding residual: CUBE_DECIMALS
# puts that near 3e-8, which is reported separately and deliberately NOT gated,
# since rounding the cube for file size is intended.
MEAN_MATCH_TOL = 1e-12

# First-principles band on 2024 GMST vs 1850-1900 (IGCC 2024 ~1.24 K over
# 2015-2024). Outside this, suspect the input, not the result.
GMST_2024_LO, GMST_2024_HI = 0.9, 2.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ssp", required=True, choices=sorted(SSP_TO_MARKER))
    ap.add_argument("--marker", default=None,
                    help="override the SSP_TO_MARKER default; 'none' selects the "
                         "unsuffixed volcanic_solar.csv (marker-free)")
    # ⚠ MARKER-FREE MODE (2026-09-02), mirroring build_fair_mean_v160.py. With
    # --calib-dir calibration_v160 --marker none, Land use runs `calculated` from the
    # scenario's own cumulative CO2 AFOLU and Irrigation is ONE shared trajectory, so
    # NO CMIP7 marker choice enters. ⚠⚠ That uses the stock posterior OFF-DESIGN: the
    # constrained parameters are byte-identical to the prescribed-land-use calibration,
    # so 1.6.0 was never constrained under `calculated`. Measured cost over the
    # constraining period is 2-6 % of the ensemble's own p5-p95 spread — smaller than the
    # marker ambiguity it removes (up to 19 %). State the deviation wherever it is used.
    ap.add_argument("--calib-dir", type=Path, required=True,
                    help="calibration directory: calibration_v160_prod for "
                         "production, calibration_v160 for marker-free")
    ap.add_argument("--emissions-dir", type=Path, default=None,
                    help="where the emissions live, if not the calib dir")
    ap.add_argument("--out-tag", default=None,
                    help="output filename infix: fair_cube_gmst_<SSP><OUT_TAG>_raw.csv")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--truncation-reference-dir", type=Path, default=None,
                    help="directory holding the quarantined pre-fix 452-row "
                         "cubes; enables the [CUBE-TRUNCATION] gate")
    args = ap.parse_args()

    ssp = args.ssp
    marker = args.marker or SSP_TO_MARKER[ssp]
    cal_dir = args.calib_dir
    emis_dir = args.emissions_dir if args.emissions_dir else cal_dir
    params = cal_dir / f"calibrated_constrained_parameters_{CALIB_VERSION}.csv"
    species_file = cal_dir / f"species_configs_properties_{CALIB_VERSION}.csv"
    forcing_file = (cal_dir / "volcanic_solar.csv") if str(marker).lower() == "none" \
        else cal_dir / f"volcanic_solar_{marker}.csv"
    emis = emis_dir / f"emissions_v160_{BASIS_SLUG}_{ssp}.csv"
    for p in (params, species_file, forcing_file, emis):
        if not p.exists():
            sys.exit(f"ERROR: missing {p}")

    provenance = (f"{MODEL_LABEL} (calib {CALIB_VERSION}); CMIP7 historical "
                  f"1750-2023 + harmonized RCMIP {ssp} future ({emis.name}); "
                  f"prescribed land use + irrigation from CMIP7 marker {marker}; "
                  f"anomaly vs {BASELINE_START}-{BASELINE_END}")

    df_configs = pd.read_csv(params, index_col=0)
    configs = df_configs.index.tolist()
    species, properties = read_properties(str(species_file))
    prescribed = [s for s in species if properties[s]["input_mode"] == "forcing"]

    print(f"── per-config cube, {ssp} ───────────────────────────────────────────")
    # ⚠ the RESOLVED dir, not the module default — otherwise a marker-free run prints
    # "calibration_v160_prod" and the provenance line lies about what it used.
    print(f"  calibration : {CALIB_VERSION}  ({cal_dir.name}, marker {marker})")
    print(f"  emissions   : {emis.name}")
    print(f"  forcing     : {forcing_file.name}  (marker {marker})")
    print(f"  {len(configs)} configs, {len(species)} species, "
          f"prescribed: {prescribed}")

    f = FAIR()
    f.define_time(START_YEAR, END_YEAR, step=1)
    f.define_scenarios([ssp])
    f.define_species(species, properties)
    f.ch4_method = "Thornhill2021"
    f.define_configs(configs)
    f.allocate()

    t = FAIR()
    t.define_time(START_YEAR, END_YEAR, step=1)
    t.define_scenarios(["all"])
    t.define_species(species, properties)
    t.ch4_method = "Thornhill2021"
    t.define_configs(configs)
    t.allocate()
    t.fill_from_csv(emissions_file=str(emis), forcing_file=str(forcing_file))
    emis_arr = t.emissions.sel(scenario="all").values
    forc_arr = t.forcing.sel(scenario="all").values
    del t

    idx = {s: i for i, s in enumerate(species)}
    blank = [s for s in prescribed if np.isnan(forc_arr[:, :, idx[s]]).all()]
    if blank:
        sys.exit(f"ERROR: {forcing_file.name} has no rows for prescribed {blank}")
    ed = [s for s in species if properties[s]["input_mode"] == "emissions"]
    miss = sorted({s for s in ed if np.isnan(emis_arr[:, :, idx[s]]).any()})
    if miss:
        sys.exit(f"ERROR: {emis.name} has NaN for emissions-driven {miss}")

    f.emissions.loc[dict(scenario=ssp)] = emis_arr
    f.forcing.loc[dict(scenario=ssp)] = forc_arr
    f.fill_species_configs(str(species_file))
    f.override_defaults(str(params))

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
    for lab, arr, mean_col, mean_name in [
            ("gmst", gmst, "gmst_C", "gmst"),
            ("ohc", ohc, "ohc_1e22J", "ohc")]:
        cube = pd.DataFrame(np.round(arr[cy, :], CUBE_DECIMALS),
                            columns=[f"cfg_{c}" for c in configs])
        cube.insert(0, "year", years[cy])
        _t = args.out_tag or ""
        cp = args.out_dir / f"fair_cube_{lab}_{ssp}{_t}_raw.csv"
        cube.to_csv(cp, index=False)
        mp = args.out_dir / f"fair_mean_{mean_name}_{ssp}{_t}.csv"
        pd.DataFrame({"year": years, mean_col: arr.mean(axis=1)}).to_csv(mp, index=False)
        # [CUBE-TRUNCATION] if a pre-fix cube is in the quarantine, this re-run
        # must reproduce it EXACTLY with its out-of-range rows dropped -- proving
        # the fix moved the year axis and nothing else. Editing the row off on
        # disk would have been quicker and would have left the artifact
        # unreproducible from its own script.
        qp = (args.truncation_reference_dir / cp.name
              if args.truncation_reference_dir else None)
        if qp is None:
            print("  [CUBE-TRUNCATION] skipped: no --truncation-reference-dir")
        elif not qp.exists():
            print(f"  [CUBE-TRUNCATION] skipped: no {cp.name} in "
                  f"{args.truncation_reference_dir}")
        else:
            old = pd.read_csv(qp).set_index("year")
            new = pd.read_csv(cp).set_index("year")
            keep = old.index[(old.index >= CUBE_FIRST_YEAR)
                             & (old.index <= CUBE_LAST_YEAR)]
            same_ax = list(new.index) == list(keep)
            w = (float(np.abs(new.loc[keep].to_numpy()
                              - old.loc[keep].to_numpy()).max())
                 if same_ax and list(new.columns) == list(old.columns) else float("inf"))
            print(f"  [CUBE-TRUNCATION] vs quarantined pre-fix: axis "
                  f"{'matches' if same_ax else 'DIFFERS'}, dropped "
                  f"{len(old.index) - len(keep)} row(s), max |diff| on the kept "
                  f"rows = {w:.3e}  {'PASS' if w == 0.0 else 'FAIL'}")
            if w != 0.0:
                sys.exit("ERROR: the re-run changed more than the year axis")
        written[lab] = (cp, mp, arr)
        print(f"  wrote {cp.name}  ({cp.stat().st_size / 1e6:.1f} MB, "
              f"{len(configs)} cfg x {int(cy.sum())} yr)")
        print(f"  wrote {mp.name}")

    # [MEAN-MATCH] the mean file IS the cube's reduction, exactly, in memory.
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
