"""
build_forcing.py -- rebuild every FaIR forcing file Ladrillo ships in
data/forcing/, from the raw inputs, in dependency order.

    FaIR 2.2.4 (calib 1.6.0); CMIP7 historical emissions 1750-2023; GMST and OHC
    as anomalies vs 1850-1900; 841 configs.

THE CHAIN (each step is one ported research builder, run as its own process):

  calib      1  setup_calibration_v160.py            Zenodo 18828694 -> calibration_v160/
             2  setup_calibration_v160_production.py -> calibration_v160_prod/
  emissions  3  build_emissions_v160_ssp245.py       CMIP7 hist + chrisroadmap SSP2-4.5 + AR6 r(t)
                                                     -> calibration_v160/emissions_v160_cmip7_ssp245_harmonized.csv
             4  build_emissions_v160_ssp_rcmip.py    CMIP7 hist + RCMIP ssp126/245/585
                                                     -> calibration_v160_prod/emissions_v160_cmip7harm_<ssp>.csv
             6  build_vanvuuren_extended_splice.py   Zenodo 20713982 v1.1.1 + vv_template_skeleton.csv
                                                     -> vanvuuren/spliced_ext_zenodo/
             7  build_emissions_v160_cmip7harm_vv.py CMIP7 hist rebase -> calibration_v160_prod/emissions_v160_cmip7harm_vv<M>.csv
  fair       8  build_fair_mean_v160.py              ssp245harm (the CALIBRATION forcing)
             9  build_fair_cube_v160.py   x3         ssp126 ssp245 ssp585 (cube + mean)
            10  build_fair_cube_vv_v160.py x7        vvVL vvLN vvL vvML vvM vvHL vvH (cube + mean)
            11  build_fair_cube_spliced_v160.py      fair_cube_gmst_ssp585_spliced.csv

There is no step 5. It built a van Vuuren template from the GATED IIASA
scenariomip-2026-prerelease file plus Smith 2024, and step 6 read that template.
Since 2026-10-08 the vv files are built on the `zenodo` variant, whose values
from 2023.5 on all come from the published Zenodo 20713982 v1.1.1 file; the
template contributed only names, units and a column grid, which
vv_template_skeleton.csv now carries. Step ids keep their numbers.

⚠ THE CALIBRATION FORCING IS A MIX OF TWO DIRECTORIES. Step 8 runs
calibration_v160_PROD (prescribed land use, marker M) on calibration_v160's
emissions file. `--calib-dir calibration_v160 --marker none` -- the obvious
reading -- is 0.030 K cooler at 2024 and reproduces neither shipped file.

Everything is written under --work-dir. data/forcing/ is touched ONLY with
--install, and only after every expected file was built in this run's plan.

Inputs come from EXACTLY ONE of
    --inputs-dir DIR        the layout fetch_inputs.py writes (see inputs.sha256)
    --fairtofredi-root DIR  the maintainer's research checkout's own copies
Every step declares the inputs it reads. Before anything runs, each input is
checked against inputs.sha256:
  * PRESENT but failing its checksum -> the build ABORTS; nothing runs;
  * ABSENT -> the steps that read it, and every step downstream of them, are
    dropped; the build says which steps and which of the 43 files will not be
    built, and builds the rest. --install then refuses.

The builders are run from a FROZEN COPY of this directory
(<work-dir>/tools_snapshot/), so editing a script mid-build cannot change a
running job. Every command, its wall time, the script and input sha256s and the
package versions go to <work-dir>/BUILD_PROVENANCE.json.

    python build_forcing.py --work-dir /tmp/forcing_build --inputs-dir INPUTS
    python check_forcing_identity.py --candidate-dir /tmp/forcing_build/forcing
    #   after a partial build (an input absent): add --allow-missing
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import forcing_inputs as fi  # noqa: E402
import setup_calibration_v160 as _cal  # noqa: E402  (its NEED list: step 1's raw files)
import setup_calibration_v160_production as _prod  # noqa: E402  (step 2's raw files)

# -- Named constants: every path, label and expected filename derives from these
MODEL_LABEL = "FaIR 2.2.4 (calib 1.6.0)"
CALIB_VERSION = "1.6.0"
LADRILLO_ROOT = HERE.parents[1]
INSTALL_DIR = LADRILLO_ROOT / "data" / "forcing"

CALIB_TAG = "ssp245harm"                 # the Ladrillo calibration forcing
CALIB_SCENARIO_LABEL = "ssp245_harmonized"
CALIB_MARKER = "M"
CALIB_EMIS_NAME = "emissions_v160_cmip7_ssp245_harmonized.csv"
SSPS = ["ssp126", "ssp245", "ssp585"]
VV_MARKERS = ["VL", "LN", "L", "ML", "M", "HL", "H"]
VV_PREFIX = "vv"
VV_TAIL_VARIANT = "zenodo"               # the variant the cmip7harm vv files rebase
VV_SKELETON_NAME = "vv_template_skeleton.csv"   # names + units for step 6 (no numbers)
SPLICED_SSP = "ssp585"
SPLICED_FIELDS = ["gmst"]                # Ladrillo ships only the GMST spliced cube
FIELDS = ["gmst", "ohc"]

STAGES = ("calib", "emissions", "fair")
SNAPSHOT_DIR = "tools_snapshot"
BLAS_ENV = {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                             "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
ZENODO_CALIB_KEYS = ("calib_params", "calib_species", "calib_volcanic",
                     "calib_solar", "calib_irrigation", "calib_landuse",
                     "cmip7_hist_emissions")
# The calib setup scripts take a DIRECTORY (--raw-dir) and read named files from
# it, so their inputs are not visible in argv. They are taken from the scripts'
# OWN file-name constants, not restated here.
STEP1_RAW_FILES = tuple(_cal.NEED)
STEP2_RAW_FILES = (_prod.VOLCANIC_NAME, _prod.SOLAR_NAME, _prod.LANDUSE_NAME,
                   _prod.IRRIGATION_NAME)
# fair.fill_from_rcmip() retrieves all three RCMIP files through pooch's cache.
RCMIP_KEYS = tuple(fi.RCMIP_FILES)


class Step(NamedTuple):
    sid: str
    stage: str
    script: str
    argv: list
    deps: list
    inputs: tuple        # forcing_inputs.INPUTS keys this step reads
    outputs: tuple       # the data/forcing/ files this step writes (finals only)


def expected_outputs():
    """The 43 files of data/forcing/, derived from the constants above."""
    out = [f"fair_mean_{f}_{CALIB_TAG}.csv" for f in FIELDS]
    for scen in SSPS + [f"{VV_PREFIX}{m}" for m in VV_MARKERS]:
        for f in FIELDS:
            out += [f"fair_cube_{f}_{scen}_raw.csv", f"fair_mean_{f}_{scen}.csv"]
    out += [f"fair_cube_{f}_{SPLICED_SSP}_spliced.csv" for f in SPLICED_FIELDS]
    return sorted(out)


def _scenario_outputs(scen):
    return tuple(n for f in FIELDS
                 for n in (f"fair_cube_{f}_{scen}_raw.csv", f"fair_mean_{f}_{scen}.csv"))


def _raw_keys(filenames):
    """Zenodo 18828694 file names -> their INPUTS keys (exactly one each)."""
    keys = []
    for name in filenames:
        hit = [k for k in ZENODO_CALIB_KEYS if Path(fi.INPUTS[k]["rel"]).name == name]
        if len(hit) != 1:
            sys.exit(f"ERROR: raw file {name!r} maps to {hit} in forcing_inputs.INPUTS")
        keys += hit
    return tuple(keys)


def plan(w, inp, args, skeleton):
    """[Step] in dependency order. Each step records the INPUTS keys it reads:
    `use(key)` resolves a path for argv AND records the key against the step
    being assembled; `add()` takes the recorded keys and resets the record.
    """
    cal, prod, vv, out = (w / "calibration_v160", w / "calibration_v160_prod",
                          w / "vanvuuren", w / "forcing")
    raw_dirs = {inp(k).parent for k in ZENODO_CALIB_KEYS}
    if len(raw_dirs) != 1:
        sys.exit(f"ERROR: the Zenodo {fi.ZENODO_CALIB_RECORD} files must share one "
                 f"directory, found {sorted(map(str, raw_dirs))}")
    raw_root = raw_dirs.pop()
    species_prod = prod / f"species_configs_properties_{CALIB_VERSION}.csv"
    s, used = [], []

    def use(key):
        used.append(key)
        return inp(key)

    def raw_dir(filenames):
        used.extend(_raw_keys(filenames))
        return raw_root

    def add(sid, stage, script, argv, deps, extra_inputs=(), outputs=()):
        keys = tuple(dict.fromkeys([*used, *extra_inputs]))
        used.clear()
        s.append(Step(sid, stage, script, [str(a) for a in argv], list(deps),
                      keys, tuple(outputs)))

    add("01_calib_v160", "calib", "setup_calibration_v160.py",
        ["--raw-dir", raw_dir(STEP1_RAW_FILES), "--out-dir", cal], [])
    add("02_calib_v160_prod", "calib", "setup_calibration_v160_production.py",
        ["--src-calib-dir", cal, "--raw-dir", raw_dir(STEP2_RAW_FILES),
         "--out-dir", prod], ["01_calib_v160"])
    add("03_emis_ssp245harm", "emissions", "build_emissions_v160_ssp245.py",
        ["--hist-file", use("cmip7_hist_emissions"), "--chris-file", use("chris_ssp245"),
         "--ar6-file", use("ar6_ssp245"), "--out-dir", cal]
        + (["--v145-builder", args.v145_builder] if args.v145_builder else []),
        ["01_calib_v160"])
    add("04_emis_ssp_rcmip", "emissions", "build_emissions_v160_ssp_rcmip.py",
        ["--ssps", *SSPS, "--hist-file", use("cmip7_hist_emissions"),
         "--species-file", species_prod, "--out-dir", prod],
        ["02_calib_v160_prod"], extra_inputs=RCMIP_KEYS)
    add("06_vv_spliced_ext", "emissions", "build_vanvuuren_extended_splice.py",
        ["--ext-csv", use("vv_ext_zenodo"), "--skeleton", skeleton,
         "--out-root", vv, "--variants", VV_TAIL_VARIANT], [])
    add("07_emis_vv_cmip7harm", "emissions", "build_emissions_v160_cmip7harm_vv.py",
        ["--hist-file", use("cmip7_hist_emissions"), "--spliced-ext-root", vv,
         "--variant", VV_TAIL_VARIANT, "--out-dir", prod],
        ["06_vv_spliced_ext", "02_calib_v160_prod"])
    add(f"08_mean_{CALIB_TAG}", "fair", "build_fair_mean_v160.py",
        ["--calib-dir", prod, "--emissions-file", cal / CALIB_EMIS_NAME,
         "--tag", CALIB_TAG, "--scenario-label", CALIB_SCENARIO_LABEL,
         "--marker", CALIB_MARKER, "--out-dir", out],
        ["02_calib_v160_prod", "03_emis_ssp245harm"],
        outputs=[f"fair_mean_{f}_{CALIB_TAG}.csv" for f in FIELDS])
    for ssp in SSPS:
        add(f"09_cube_{ssp}", "fair", "build_fair_cube_v160.py",
            ["--ssp", ssp, "--calib-dir", prod, "--out-dir", out]
            + (["--truncation-reference-dir", args.truncation_reference_dir]
               if args.truncation_reference_dir else []),
            ["04_emis_ssp_rcmip"], outputs=_scenario_outputs(ssp))
    for m in VV_MARKERS:
        add(f"10_cube_{VV_PREFIX}{m}", "fair", "build_fair_cube_vv_v160.py",
            ["--marker", m, "--calib-dir", prod, "--out-dir", out],
            ["07_emis_vv_cmip7harm"], outputs=_scenario_outputs(f"{VV_PREFIX}{m}"))
    add(f"11_spliced_{SPLICED_SSP}", "fair", "build_fair_cube_spliced_v160.py",
        ["--ssp", SPLICED_SSP, "--fields", ",".join(SPLICED_FIELDS),
         "--in-dir", out, "--out-dir", out]
        + (["--port-gate-dir", args.spliced_port_gate_dir]
           if args.spliced_port_gate_dir else ["--skip-port-gate"]),
        [f"09_cube_{SPLICED_SSP}"],
        outputs=[f"fair_cube_{f}_{SPLICED_SSP}_spliced.csv" for f in SPLICED_FIELDS])
    _check_plan(s)
    return s


def _check_plan(steps):
    """The plan, the input table and the expected file list must agree."""
    ids = [x.sid for x in steps]
    bad_dep = [(x.sid, d) for x in steps for d in x.deps if d not in ids[:ids.index(x.sid)]]
    if bad_dep:
        sys.exit(f"ERROR: steps depend on unknown or later steps: {bad_dep}")
    outs = [n for x in steps for n in x.outputs]
    if sorted(outs) != expected_outputs():
        sys.exit("ERROR: the steps' declared outputs are not the 43 expected files: "
                 f"extra {sorted(set(outs) - set(expected_outputs()))}, missing "
                 f"{sorted(set(expected_outputs()) - set(outs))}")
    unused = sorted(set(fi.INPUTS) - {k for x in steps for k in x.inputs})
    if unused:
        sys.exit(f"ERROR: forcing_inputs.INPUTS lists inputs no step reads: {unused}")


def prune(steps, absent):
    """-> (kept steps, {dropped sid: reason}). A step is dropped if it reads an
    absent input or depends on a dropped step."""
    dropped = {}
    for x in steps:                                    # steps are in dependency order
        miss = [k for k in x.inputs if k in absent]
        up = [d for d in x.deps if d in dropped]
        if miss:
            dropped[x.sid] = f"input absent: {', '.join(miss)}"
        elif up:
            dropped[x.sid] = f"depends on dropped {', '.join(up)}"
    return [x for x in steps if x.sid not in dropped], dropped


def run_step(step, snap, logs, python):
    cmd = [python, str(snap / step.script), *step.argv]
    env = {**os.environ, **BLAS_ENV}
    t0 = time.time()
    with open(logs / f"{step.sid}.log", "w") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        log.flush()
        rc = subprocess.run(cmd, cwd=snap, env=env, stdout=log,
                            stderr=subprocess.STDOUT).returncode
    return step.sid, rc, time.time() - t0, cmd


def freeze_tools(dest):
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    hashes = {}
    for p in sorted(HERE.glob("*.py")) + [HERE / VV_SKELETON_NAME]:
        shutil.copy2(p, dest / p.name)
        hashes[p.name] = fi.sha256(p)
    shutil.copy2(fi.SHA_FILE, dest / fi.SHA_FILE.name)
    return hashes


def env_versions():
    import importlib.metadata as md
    pk = {}
    for name in ("fair", "numpy", "pandas", "xarray", "scipy", "pooch", "tqdm"):
        try:
            pk[name] = md.version(name)
        except md.PackageNotFoundError:
            pk[name] = None
    return {"python": sys.version.split()[0], "platform": platform.platform(),
            "packages": pk}


def git_state():
    try:
        head = subprocess.run(["git", "-C", str(LADRILLO_ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(LADRILLO_ROOT), "status",
                                "--porcelain", "--", str(HERE)],
                               capture_output=True, text=True).stdout.strip()
        return {"ladrillo_head": head, "tools_forcing_dirty": bool(dirty)}
    except OSError:
        return {}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work-dir", type=Path, required=True,
                    help="scratch directory; finals land in <work-dir>/forcing/")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--inputs-dir", type=Path,
                     help="inputs laid out as fetch_inputs.py writes them")
    src.add_argument("--fairtofredi-root", type=Path,
                     help="read the maintainer's research checkout's own copies instead")
    ap.add_argument("--stages", nargs="+", default=list(STAGES), choices=STAGES)
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel FaIR runs (each single-threaded; ~1-2 GB each)")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--install", action="store_true",
                    help=f"copy the {len(expected_outputs())} finished files into "
                         f"{INSTALL_DIR} (OVERWRITES the shipped forcing)")
    ap.add_argument("--v145-builder", type=Path, default=None,
                    help="research build_emissions_v145_ssp245.py: cross-check "
                         "the vendored species maps")
    ap.add_argument("--truncation-reference-dir", type=Path, default=None,
                    help="research quarantine 20260828_cube_year_2301 (optional gate)")
    ap.add_argument("--spliced-port-gate-dir", type=Path, default=None,
                    help="research quarantine 20260828_calib160_migration "
                         "(optional gate; skipped and SAID so without it)")
    args = ap.parse_args()

    w = args.work_dir.resolve()
    w.mkdir(parents=True, exist_ok=True)
    logs = w / "logs"
    logs.mkdir(exist_ok=True)
    snap = w / SNAPSHOT_DIR

    def inp(key):
        return fi.resolve(key, args.inputs_dir, args.fairtofredi_root)

    print(f"build_forcing.py -- {MODEL_LABEL}")
    print(f"  work dir : {w}")
    print(f"  inputs   : {args.inputs_dir or args.fairtofredi_root} "
          f"({'--inputs-dir' if args.inputs_dir else '--fairtofredi-root'})")
    print(f"  stages   : {' '.join(args.stages)}   jobs: {args.jobs}")

    full = plan(w, inp, args, snap / VV_SKELETON_NAME)
    print("\n[INPUTS] sha256 against inputs.sha256 (RCMIP also against fair's md5):")
    status = fi.check_inputs(args.inputs_dir, args.fairtofredi_root)
    bad = {k: why for k, (st, why) in status.items() if st == "bad"}
    absent = {k: why for k, (st, why) in status.items() if st == "absent"}
    if bad:
        sys.exit(f"ERROR: {len(bad)} input(s) PRESENT but failing their checksum: "
                 f"{sorted(bad)}; nothing was run.")

    kept, dropped = prune(full, absent)
    not_built = sorted(n for x in full if x.sid in dropped for n in x.outputs)
    if dropped:
        print(f"\n[PARTIAL] {len(absent)} input(s) ABSENT: {', '.join(sorted(absent))}")
        if any(k in fi.RCMIP_FILES for k in absent):
            print("  (RCMIP is read from pooch's cache; fetch_inputs.py --download primes it)")
        print(f"  {len(dropped)} step(s) will NOT run:")
        for sid, why in dropped.items():
            print(f"    {sid:<24} {why}")
        print(f"  {len(not_built)} of {len(expected_outputs())} expected files will "
              f"NOT be built:")
        for n in not_built:
            print(f"    {n}")
        print("  Everything else is built. --install will refuse.")

    steps = [x for x in kept if x.stage in args.stages]
    selected = {x.sid for x in steps}

    script_sha = freeze_tools(snap)
    print(f"\n[FROZEN] {len(script_sha)} tool files copied to {snap}")

    prov = {"model": MODEL_LABEL, "started_utc": datetime.now(timezone.utc).isoformat(),
            "work_dir": str(w), "inputs_root": str(args.inputs_dir or args.fairtofredi_root),
            "inputs_mode": "inputs-dir" if args.inputs_dir else "fairtofredi-root",
            "input_sha256": {k: fi.sha256(inp(k)) for k in fi.INPUTS if k not in absent},
            "inputs_absent": absent, "steps_dropped": dropped,
            "outputs_not_built": not_built,
            "script_sha256": script_sha, "env": env_versions(), "git": git_state(),
            "blas_env": BLAS_ENV, "steps": []}

    done, t_all = {}, time.time()
    serial = [x for x in steps if x.stage != "fair"]
    fair = [x for x in steps if x.stage == "fair"]
    for step in serial:
        print(f"\n[{step.sid}] {step.script} ...", flush=True)
        sid, rc, dt, cmd = run_step(step, snap, logs, args.python)
        prov["steps"].append({"id": sid, "rc": rc, "seconds": round(dt, 1), "cmd": cmd})
        print(f"  rc={rc}  {dt:.1f} s  (log: logs/{sid}.log)")
        if rc != 0:
            _write_prov(w, prov)
            sys.exit(f"ERROR: {sid} failed; see {logs / (sid + '.log')}")
        done[sid] = True

    pending = list(fair)
    with cf.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        running = {}
        while pending or running:
            ready = [x for x in pending
                     if all(d in done or d not in selected for d in x.deps)]
            for x in ready:
                if len(running) >= max(1, args.jobs):
                    break
                pending.remove(x)
                print(f"\n[{x.sid}] {x.script} {' '.join(x.argv[:2])} ... (started)",
                      flush=True)
                running[ex.submit(run_step, x, snap, logs, args.python)] = x
            if not running:
                sys.exit(f"ERROR: unsatisfiable dependencies: {[x.sid for x in pending]}")
            fin, _ = cf.wait(list(running), return_when=cf.FIRST_COMPLETED)
            for fut in fin:
                running.pop(fut)
                sid, rc, dt, cmd = fut.result()
                prov["steps"].append({"id": sid, "rc": rc, "seconds": round(dt, 1),
                                      "cmd": cmd})
                print(f"[{sid}] rc={rc}  {dt:.1f} s", flush=True)
                if rc != 0:
                    _write_prov(w, prov)
                    for f2 in running:
                        f2.cancel()
                    sys.exit(f"ERROR: {sid} failed; see {logs / (sid + '.log')}")
                done[sid] = True

    prov["wall_seconds"] = round(time.time() - t_all, 1)
    out = w / "forcing"
    want = expected_outputs()
    have = sorted(p.name for p in out.glob("*.csv")) if out.exists() else []
    planned = sorted(n for x in steps for n in x.outputs)
    failed = sorted(set(planned) - set(have))      # a step ran and did not write it
    missing = sorted(set(want) - set(have))
    extra = sorted(set(have) - set(want))
    prov["outputs"] = {n: fi.sha256(out / n) for n in want if (out / n).exists()}
    prov["missing_outputs"], prov["extra_outputs"] = missing, extra
    _write_prov(w, prov)
    print(f"\nDone in {prov['wall_seconds']:.0f} s. {len(want) - len(missing)} of "
          f"{len(want)} expected files in {out}"
          + (f"; extra (not shipped) {extra}" if extra else ""))
    if not_built:
        print(f"PARTIAL BUILD: {len(not_built)} file(s) not built because an input "
              f"is absent: {not_built}")
    print(f"Provenance: {w / 'BUILD_PROVENANCE.json'}")
    print(f"Next: python {HERE / 'check_forcing_identity.py'} --candidate-dir {out}"
          + (" --allow-missing" if missing else ""))
    if failed:
        sys.exit(f"ERROR: steps ran but these expected files are missing: {failed}")

    if args.install:
        if dropped or missing:
            sys.exit("ERROR: --install refused: not every expected file was built "
                     f"({len(dropped)} step(s) dropped, {len(missing)} file(s) missing).")
        for n in want:
            shutil.copy2(out / n, INSTALL_DIR / n)
        print(f"\n[INSTALL] copied {len(want)} files into {INSTALL_DIR}. "
              f"Re-run check_forcing_identity.py and update data/MANIFEST.csv.")


def _write_prov(w, prov):
    (w / "BUILD_PROVENANCE.json").write_text(json.dumps(prov, indent=1, default=str))


if __name__ == "__main__":
    main()
