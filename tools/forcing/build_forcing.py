"""
build_forcing.py -- rebuild every FaIR forcing file Ladrillo ships in
data/forcing/, from the raw inputs, in dependency order.

    FaIR 2.2.4 (calib 1.6.0); CMIP7 historical emissions 1750-2023; GMST and OHC
    as anomalies vs 1850-1900; 841 configs.

THE CHAIN (each step is one ported research builder, run as its own process):

  calib      1  setup_calibration_v160.py            Zenodo 18828694 -> calibration_v160/
             2  setup_calibration_v160_production.py -> calibration_v160_prod/
  emissions  3  build_emissions_v160_ssp245.py       CMIP7 hist + chrisroadmap SSP2-4.5
                                                     -> calibration_v160/emissions_v160_cmip7_ssp245_harmonized.csv
             4  build_emissions_v160_ssp_rcmip.py    CMIP7 hist + RCMIP ssp126/245/585
                                                     -> calibration_v160_prod/emissions_v160_cmip7harm_<ssp>.csv
             5  build_vanvuuren_spliced.py           Smith 2024 + IIASA prerelease -> vanvuuren/spliced/
             6  build_vanvuuren_extended_splice.py   + Zenodo 20713982 v1.1.1 tail -> vanvuuren/spliced_ext_harmonized/
             7  build_emissions_v160_cmip7harm_vv.py CMIP7 hist rebase -> calibration_v160_prod/emissions_v160_cmip7harm_vv<M>.csv
  fair       8  build_fair_mean_v160.py              ssp245harm (the CALIBRATION forcing)
             9  build_fair_cube_v160.py   x3         ssp126 ssp245 ssp585 (cube + mean)
            10  build_fair_cube_vv_v160.py x7        vvVL vvLN vvL vvML vvM vvHL vvH (cube + mean)
            11  build_fair_cube_spliced_v160.py      fair_cube_gmst_ssp585_spliced.csv

⚠ THE CALIBRATION FORCING IS A MIX OF TWO DIRECTORIES. Step 8 runs
calibration_v160_PROD (prescribed land use, marker M) on calibration_v160's
emissions file. `--calib-dir calibration_v160 --marker none` -- the obvious
reading -- is 0.030 K cooler at 2024 and reproduces neither shipped file.

Everything is written under --work-dir. data/forcing/ is touched ONLY with
--install, and only after every expected file exists.

Inputs come from EXACTLY ONE of
    --inputs-dir DIR        the layout fetch_inputs.py writes (see inputs.sha256)
    --fairtofredi-root DIR  the research repo's own copies (Marcus's machine)
and every input is checked against inputs.sha256 before anything runs.

The builders are run from a FROZEN COPY of this directory
(<work-dir>/tools_snapshot/), so editing a script mid-build cannot change a
running job. Every command, its wall time, the script and input sha256s and the
package versions go to <work-dir>/BUILD_PROVENANCE.json.

    python build_forcing.py --work-dir /tmp/forcing_build --inputs-dir INPUTS
    python check_forcing_identity.py --candidate-dir /tmp/forcing_build/forcing
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import forcing_inputs as fi  # noqa: E402

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
VV_TAIL_VARIANT = "harmonized"           # the variant the cmip7harm vv files rebase
SPLICED_SSP = "ssp585"
SPLICED_FIELDS = ["gmst"]                # Ladrillo ships only the GMST spliced cube
FIELDS = ["gmst", "ohc"]

STAGES = ("calib", "emissions", "fair")
BLAS_ENV = {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                             "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
ZENODO_CALIB_KEYS = ("calib_params", "calib_species", "calib_volcanic",
                     "calib_solar", "calib_irrigation", "calib_landuse",
                     "cmip7_hist_emissions")


def expected_outputs():
    """The 43 files of data/forcing/, derived from the constants above."""
    out = [f"fair_mean_{f}_{CALIB_TAG}.csv" for f in FIELDS]
    for scen in SSPS + [f"{VV_PREFIX}{m}" for m in VV_MARKERS]:
        for f in FIELDS:
            out += [f"fair_cube_{f}_{scen}_raw.csv", f"fair_mean_{f}_{scen}.csv"]
    out += [f"fair_cube_{f}_{SPLICED_SSP}_spliced.csv" for f in SPLICED_FIELDS]
    return sorted(out)


def plan(w, inp, args):
    """[(step_id, stage, script, argv, depends_on)] in dependency order."""
    cal, prod, vv, out = (w / "calibration_v160", w / "calibration_v160_prod",
                          w / "vanvuuren", w / "forcing")
    raw = {k: inp(k) for k in ZENODO_CALIB_KEYS}
    raw_dirs = {p.parent for p in raw.values()}
    if len(raw_dirs) != 1:
        sys.exit(f"ERROR: the Zenodo {fi.ZENODO_CALIB_RECORD} files must share one "
                 f"directory, found {sorted(map(str, raw_dirs))}")
    raw_dir = raw_dirs.pop()
    hist = raw["cmip7_hist_emissions"]
    species_prod = prod / f"species_configs_properties_{CALIB_VERSION}.csv"
    s = []
    s.append(("01_calib_v160", "calib", "setup_calibration_v160.py",
              ["--raw-dir", raw_dir, "--out-dir", cal], []))
    s.append(("02_calib_v160_prod", "calib", "setup_calibration_v160_production.py",
              ["--src-calib-dir", cal, "--raw-dir", raw_dir, "--out-dir", prod],
              ["01_calib_v160"]))
    s.append(("03_emis_ssp245harm", "emissions", "build_emissions_v160_ssp245.py",
              ["--hist-file", hist, "--chris-file", inp("chris_ssp245"),
               "--ar6-file", inp("ar6_ssp245"), "--out-dir", cal]
              + (["--v145-builder", args.v145_builder] if args.v145_builder else []),
              ["01_calib_v160"]))
    s.append(("04_emis_ssp_rcmip", "emissions", "build_emissions_v160_ssp_rcmip.py",
              ["--ssps", *SSPS, "--hist-file", hist, "--species-file", species_prod,
               "--out-dir", prod], ["02_calib_v160_prod"]))
    s.append(("05_vv_spliced", "emissions", "build_vanvuuren_spliced.py",
              ["--vv-file", inp("vv_iiasa_prerelease"),
               "--smith-file", inp("smith2024_history"), "--out-dir", vv / "spliced"],
              []))
    s.append(("06_vv_spliced_ext", "emissions", "build_vanvuuren_extended_splice.py",
              ["--ext-csv", inp("vv_ext_zenodo"), "--spliced-dir", vv / "spliced",
               "--out-root", vv, "--variants", VV_TAIL_VARIANT], ["05_vv_spliced"]))
    s.append(("07_emis_vv_cmip7harm", "emissions", "build_emissions_v160_cmip7harm_vv.py",
              ["--hist-file", hist, "--spliced-ext-root", vv, "--variant",
               VV_TAIL_VARIANT, "--out-dir", prod],
              ["06_vv_spliced_ext", "02_calib_v160_prod"]))
    s.append((f"08_mean_{CALIB_TAG}", "fair", "build_fair_mean_v160.py",
              ["--calib-dir", prod, "--emissions-file", cal / CALIB_EMIS_NAME,
               "--tag", CALIB_TAG, "--scenario-label", CALIB_SCENARIO_LABEL,
               "--marker", CALIB_MARKER, "--out-dir", out],
              ["02_calib_v160_prod", "03_emis_ssp245harm"]))
    for ssp in SSPS:
        s.append((f"09_cube_{ssp}", "fair", "build_fair_cube_v160.py",
                  ["--ssp", ssp, "--calib-dir", prod, "--out-dir", out]
                  + (["--truncation-reference-dir", args.truncation_reference_dir]
                     if args.truncation_reference_dir else []),
                  ["04_emis_ssp_rcmip"]))
    for m in VV_MARKERS:
        s.append((f"10_cube_{VV_PREFIX}{m}", "fair", "build_fair_cube_vv_v160.py",
                  ["--marker", m, "--calib-dir", prod, "--out-dir", out],
                  ["07_emis_vv_cmip7harm"]))
    s.append((f"11_spliced_{SPLICED_SSP}", "fair", "build_fair_cube_spliced_v160.py",
              ["--ssp", SPLICED_SSP, "--fields", ",".join(SPLICED_FIELDS),
               "--in-dir", out, "--out-dir", out]
              + (["--port-gate-dir", args.spliced_port_gate_dir]
                 if args.spliced_port_gate_dir else ["--skip-port-gate"]),
              [f"09_cube_{SPLICED_SSP}"]))
    return [(sid, st, sc, [str(a) for a in av], dep) for sid, st, sc, av, dep in s]


def run_step(step, snap, logs, python):
    sid, _stage, script, argv, _dep = step
    cmd = [python, str(snap / script), *argv]
    env = {**os.environ, **BLAS_ENV}
    t0 = time.time()
    with open(logs / f"{sid}.log", "w") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        log.flush()
        rc = subprocess.run(cmd, cwd=snap, env=env, stdout=log,
                            stderr=subprocess.STDOUT).returncode
    return sid, rc, time.time() - t0, cmd


def freeze_tools(dest):
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    hashes = {}
    for p in sorted(HERE.glob("*.py")):
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
                     help="read the research repo's own copies instead")
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

    def inp(key):
        return fi.resolve(key, args.inputs_dir, args.fairtofredi_root)

    print(f"build_forcing.py -- {MODEL_LABEL}")
    print(f"  work dir : {w}")
    print(f"  inputs   : {args.inputs_dir or args.fairtofredi_root} "
          f"({'--inputs-dir' if args.inputs_dir else '--fairtofredi-root'})")
    print(f"  stages   : {' '.join(args.stages)}   jobs: {args.jobs}")
    print("\n[INPUTS] sha256 against inputs.sha256:")
    fails = fi.verify(args.inputs_dir, args.fairtofredi_root)
    if fails:
        for k, why in fails:
            print(f"  FAIL {k}: {why}")
        sys.exit("ERROR: inputs do not match inputs.sha256; nothing was run.")
    fi.require_rcmip_cache()

    snap = w / "tools_snapshot"
    script_sha = freeze_tools(snap)
    print(f"\n[FROZEN] {len(script_sha)} scripts copied to {snap}")

    steps = [s for s in plan(w, inp, args) if s[1] in args.stages]
    prov = {"model": MODEL_LABEL, "started_utc": datetime.now(timezone.utc).isoformat(),
            "work_dir": str(w), "inputs_root": str(args.inputs_dir or args.fairtofredi_root),
            "inputs_mode": "inputs-dir" if args.inputs_dir else "fairtofredi-root",
            "input_sha256": {k: fi.sha256(inp(k)) for k in fi.INPUTS},
            "script_sha256": script_sha, "env": env_versions(), "git": git_state(),
            "blas_env": BLAS_ENV, "steps": []}

    done, t_all = {}, time.time()
    serial = [s for s in steps if s[1] != "fair"]
    fair = [s for s in steps if s[1] == "fair"]
    for step in serial:
        print(f"\n[{step[0]}] {step[2]} ...", flush=True)
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
            ready = [s for s in pending
                     if all(d in done or not any(d == x[0] for x in steps)
                            for d in s[4])]
            for s in ready:
                if len(running) >= max(1, args.jobs):
                    break
                pending.remove(s)
                print(f"\n[{s[0]}] {s[2]} {' '.join(s[3][:2])} ... (started)", flush=True)
                running[ex.submit(run_step, s, snap, logs, args.python)] = s
            if not running:
                sys.exit(f"ERROR: unsatisfiable dependencies: {[s[0] for s in pending]}")
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
    missing = sorted(set(want) - set(have))
    extra = sorted(set(have) - set(want))
    prov["outputs"] = {n: fi.sha256(out / n) for n in want if (out / n).exists()}
    prov["missing_outputs"], prov["extra_outputs"] = missing, extra
    _write_prov(w, prov)
    print(f"\nDone in {prov['wall_seconds']:.0f} s. {len(want) - len(missing)} of "
          f"{len(want)} expected files in {out}"
          + (f"; MISSING {missing}" if missing and "fair" in args.stages else "")
          + (f"; extra (not shipped) {extra}" if extra else ""))
    print(f"Provenance: {w / 'BUILD_PROVENANCE.json'}")
    print(f"Next: python {HERE / 'check_forcing_identity.py'} --candidate-dir {out}")

    if args.install:
        if missing:
            sys.exit("ERROR: --install refused: not every expected file was built.")
        for n in want:
            shutil.copy2(out / n, INSTALL_DIR / n)
        print(f"\n[INSTALL] copied {len(want)} files into {INSTALL_DIR}. "
              f"Re-run check_forcing_identity.py and update data/MANIFEST.csv.")


def _write_prov(w, prov):
    (w / "BUILD_PROVENANCE.json").write_text(json.dumps(prov, indent=1, default=str))


if __name__ == "__main__":
    main()
