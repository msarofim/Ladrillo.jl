"""
build_fair_cube_spliced_v160.py

Emit the SPLICED companion to a raw per-config cube:
    fair_cube_{gmst,ohc}_<SSP>_spliced.csv

WHAT THIS IS FOR. `julia/scope_slr_fair_uncertainty.jl` reads the RAW cube and
applies the splice convention itself, in Julia, so it can be switched and gated.
The `_spliced` file is its **[SPLICE-MATCH] reference**: a Python rendering of the
same convention, which the Julia splice must reproduce (`SPLICE_MATCH_TOL = 2e-6`,
a hard `@assert`). The gate exists to prove "raw cube + convention in Julia" is a
refactor of the Python path and not a silent change of input.

WHY IT HAD TO BE REGENERATED WITH THE MIGRATION. That gate compares two things
that must share a vintage. After the calib 1.6.0 + CMIP7 migration moved the raw
cubes, a stale 1.4.5 reference made it compare a Julia splice of v160 data
against a Python splice of v145 data -- an ssp585 gap of ~0.38 K at 2100, some
190,000x the tolerance. The gate would have failed correctly, for the wrong
reason, and taken the ssp585 joint band with it. Flagged by Ladrillo's own
gate, which was not weakened -- correctly.

WHY THIS IS *NOT* PART OF THE COULON ARM, WHICH IS STILL STALE ON 1.4.5. In
`run_fair_ssp_spread.py` the splice is computed for ALL configs first (line 242)
and the Coulon arm SELECTS from it afterwards (line 246, `warm_sel =
spliced[isel, :] - ours_ref`, against Coulon's +12.0/+14.5/+17.0 degC at the
posterior-median `ais_gmst_amp`). Splice and selection are separable: the splice
is mechanical and needs only the raw cube plus the ensemble-mean path, while the
selection needs a moved FaIR ensemble AND a moved L14. This script does only the
former. `fair_coulon_arm_*`, `fair_pctile_gmst_*` and the hot-tail tables remain
stale and still need the Coulon re-run.

NO FaIR RUN. The convention is a pure function of files this migration already
produced, so this reads them off disk rather than re-integrating.

THE PORT GATE. Before writing anything, this reproduces the OLD spliced reference
from the OLD raw cube and OLD mean (both in the migration quarantine) and
requires agreement to the same 2e-6 the Julia gate uses. That proves this file's
reading of the convention matches the one that produced the shipped reference,
BEFORE it is trusted to generate the replacement. Without it, a
misremembered convention would be written straight into the reference the gate
trusts, and the gate would then certify the error.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_fair_cube_spliced_v160.py @ 8984a13. `splice()` and the write
are copied UNCHANGED; only the I/O plumbing differs:
  * the raw cube + mean are read from --in-dir and the spliced cube is written
    to --out-dir (both were SLR-RFF-BRICK/data/observations).
  * the PORT GATE's v145 reference triple (a research-repo quarantine) is read
    from --port-gate-dir. The original REFUSED to write without it; the port
    keeps that refusal unless --skip-port-gate is passed explicitly, in which
    case the skip is printed. A reviewer outside the research repo has no v145
    triple, so build_forcing.py passes --skip-port-gate unless it is given one.
Normally driven by build_forcing.py. By hand:
    python build_fair_cube_spliced_v160.py --ssp ssp585 --fields gmst \
        --in-dir OUT --out-dir OUT --skip-port-gate
"""

from pathlib import Path
import argparse
import sys

import numpy as np
import pandas as pd

# ── Named constants: these MUST match julia/scope_slr_fair_uncertainty.jl ────
SPLICE_YEAR = 2014            # last year taken from the ensemble-mean path
REF_START, REF_END = 1995, 2014   # the per-config re-referencing window (Coulon's)
CUBE_DECIMALS = 6             # both cubes are written at 6 dp

# Tolerance for the port gate, taken verbatim from the Julia gate it must match.
# Derived there from the quantity it tests: three 6-dp roundings at half-ulp
# 5e-7 each, worst case 1.5e-6, gate just above at 2e-6. A real convention
# change shows at 1e-2, so it still discriminates by ~4 orders of magnitude.
SPLICE_MATCH_TOL = 2e-6

# The research repo's v145 reference triple lives in
# SLR-RFF-BRICK/outputs/quarantine/20260828_calib160_migration; pass it as
# --port-gate-dir.

FIELDS = {"gmst": "gmst_C", "ohc": "ohc_1e22J"}


def splice(raw_cube, mean_series):
    """The convention, verbatim: mean path to SPLICE_YEAR, then the config's own
    anomaly re-referenced to REF_START-REF_END and added to the mean's."""
    years = raw_cube.index.to_numpy()
    ref = (years >= REF_START) & (years <= REF_END)
    m = mean_series.reindex(years).to_numpy(float)
    if np.isnan(m).any():
        sys.exit("ERROR: the mean file does not cover the cube's year range")
    raw = raw_cube.to_numpy(float)
    mref = m[ref].mean()
    cref = raw[ref, :].mean(axis=0)
    out = np.where((years <= SPLICE_YEAR)[:, None], m[:, None],
                   mref + (raw - cref[None, :]))
    return pd.DataFrame(out, index=raw_cube.index, columns=raw_cube.columns)


def load(cube_path, mean_path, mean_col):
    cube = pd.read_csv(cube_path).set_index("year")
    mean = pd.read_csv(mean_path).set_index("year")[mean_col]
    return cube, mean


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ssp", default="ssp585")
    ap.add_argument("--fields", default="gmst,ohc")
    ap.add_argument("--in-dir", type=Path, required=True,
                    help="directory holding fair_cube_<f>_<SSP>_raw.csv and "
                         "fair_mean_<f>_<SSP>.csv")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--port-gate-dir", type=Path, default=None,
                    help="directory holding the v145 reference triple for the "
                         "PORT GATE")
    ap.add_argument("--skip-port-gate", action="store_true",
                    help="write without the PORT GATE (no v145 triple available)")
    args = ap.parse_args()
    OBS, QUARANTINE = args.in_dir, args.port_gate_dir
    ssp = args.ssp
    fields = [f.strip() for f in args.fields.split(",")]

    print("── spliced cube (the [SPLICE-MATCH] reference) ──────────────────────")
    print(f"  ssp {ssp}; splice at {SPLICE_YEAR}; re-reference "
          f"{REF_START}-{REF_END}; {CUBE_DECIMALS} dp")

    # ── PORT GATE: reproduce the OLD reference from the OLD inputs ───────────
    print("\n  PORT GATE -- reproduce the shipped v145 reference from v145 inputs:")
    gated = False
    if QUARANTINE is None:
        if not args.skip_port_gate:
            sys.exit("ERROR: no --port-gate-dir; refusing to write a reference "
                     "this script has not been verified to produce. Pass "
                     "--skip-port-gate to write it anyway.")
        print("    SKIPPED by --skip-port-gate (no v145 reference triple supplied)")
        fields_to_gate, gated = [], True
    else:
        fields_to_gate = fields
    for fld in fields_to_gate:
        col = FIELDS[fld]
        # The v145 reference lives in the QUARANTINE, not at the canonical path.
        # It was read from OBS originally, which worked exactly once: the first
        # run overwrote its own reference with the v160 replacement, so a second
        # run gated a v145-derived splice against a v160 file. That is caught
        # (the config sets are disjoint) rather than silent, but it is still the
        # wrong source -- a gate must not read a path its own script writes.
        old_ref = QUARANTINE / f"fair_cube_{fld}_{ssp}_spliced.csv"
        old_raw = QUARANTINE / f"fair_cube_{fld}_{ssp}_raw.csv"
        old_mean = QUARANTINE / f"fair_mean_{fld}_{ssp}.csv"
        if not (old_ref.exists() and old_raw.exists() and old_mean.exists()):
            print(f"    {fld}: no v145 reference triple on disk -- CANNOT GATE")
            continue
        cube, mean = load(old_raw, old_mean, col)
        got = splice(cube, mean)
        ref = pd.read_csv(old_ref).set_index("year")
        common = [c for c in got.columns if c in ref.columns]
        if len(common) != len(ref.columns):
            sys.exit(f"ERROR: {fld} config sets differ between the v145 raw cube "
                     f"and its spliced reference")
        yy = ref.index.intersection(got.index)
        w = float(np.max(np.abs(got.loc[yy, common].to_numpy()
                                - ref.loc[yy, common].to_numpy())))
        ok = w <= SPLICE_MATCH_TOL
        print(f"    {fld}: max |this splice - shipped v145 reference| = {w:.3e} "
              f"(tol {SPLICE_MATCH_TOL:.0e})  {'PASS' if ok else 'FAIL'}")
        if not ok:
            sys.exit(f"ERROR: this script's reading of the splice convention does "
                     f"not reproduce the shipped v145 {fld} reference. Do NOT "
                     f"write a new reference from it.")
        gated = True
    if not gated:
        sys.exit("ERROR: the port gate could not run for any field; refusing to "
                 "write a reference this script has not been verified to produce.")

    # ── produce the v160 replacement ────────────────────────────────────────
    print("\n  writing the v160 replacement:")
    for fld in fields:
        col = FIELDS[fld]
        raw_p = OBS / f"fair_cube_{fld}_{ssp}_raw.csv"
        mean_p = OBS / f"fair_mean_{fld}_{ssp}.csv"
        for p in (raw_p, mean_p):
            if not p.exists():
                sys.exit(f"ERROR: missing {p}")
        cube, mean = load(raw_p, mean_p, col)
        out = splice(cube, mean).round(CUBE_DECIMALS)
        out.insert(0, "year", out.index)
        args.out_dir.mkdir(parents=True, exist_ok=True)
        dest = args.out_dir / f"fair_cube_{fld}_{ssp}_spliced.csv"
        out.to_csv(dest, index=False)

        # [SPLICE-INERT] the hindcast must be the mean path EXACTLY, pre-splice.
        # That is the entire claim the controlled experiment rests on.
        chk = pd.read_csv(dest).set_index("year")
        pre = chk.index <= SPLICE_YEAR
        m = mean.reindex(chk.index).to_numpy(float)
        worst = float(np.max(np.abs(chk.to_numpy()[pre, :]
                                    - np.round(m[pre, None], CUBE_DECIMALS))))
        print(f"    {fld}: wrote {dest.name} ({dest.stat().st_size / 1e6:.1f} MB); "
              f"[SPLICE-INERT] max |pre-{SPLICE_YEAR} - mean path| = {worst:.3e}  "
              f"{'PASS' if worst == 0.0 else 'FAIL'}")
        if worst != 0.0:
            sys.exit("ERROR: the spliced cube moves the hindcast")

    print(f"\nPROVENANCE: FaIR 2.2.4 (calib 1.6.0); spliced from "
          f"fair_cube_*_{ssp}_raw.csv + fair_mean_*_{ssp}.csv, both regenerated "
          f"on CMIP7 history; splice {SPLICE_YEAR}, re-reference "
          f"{REF_START}-{REF_END}. NOT a Coulon-arm product.")


if __name__ == "__main__":
    main()
