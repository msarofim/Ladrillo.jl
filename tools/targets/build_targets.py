#!/usr/bin/env python3
"""
build_targets.py -- rebuild Ladrillo's calibration targets and drivers from
their sources, in one command.

    python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs ...   # see its docstring
    python tools/targets/build_targets.py --inputs-dir ~/ladrillo_inputs --work-dir /tmp/targets_build
    python tools/targets/check_targets_identity.py --work-dir /tmp/targets_build

Builds, under --work-dir, the eight files Ladrillo ships that come from
observations (data/observations/calibration_targets.csv and _gsic_adj.csv,
data/drivers/glacier_block_temperature.csv and greenland_zone_temperature.csv,
data/constants/glacier_block_constants.csv, gis_amp_shape.csv,
gis_amp_shape_meta.csv, data/calibration/gis_amp_prior.csv), plus the four
intermediates the research repository also keeps. It reads data/sources/
(shipped) and the files Ladrillo does not ship from --inputs-dir (laid out by
fetch_inputs.py), merged into <work-dir>/sources/. With --from-raw it first
re-runs every maintainer extraction (tools/targets/extract/) from the raw
gridded products in <inputs-dir>/raw/ and builds from those instead.

Before any step runs, every input is checked against tools/targets/inputs.sha256;
one mismatch or missing file aborts the build, and the message says which. Each step runs from a frozen copy
of tools/targets/ (<work-dir>/tools_snapshot/), single-threaded BLAS, and the
build writes <work-dir>/BUILD_PROVENANCE.json: the commands, wall times, the
sha256 of every script, input and output, package versions, the Ladrillo
commit, and the one seed the chain uses.

It never touches data/ unless --install is given (and then only after every
file was built). Install only on purpose: the shipped posterior was calibrated
on the shipped targets.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import target_inputs as TI  # noqa: E402

REPO = TI.REPO_ROOT
BLAS_ENV = {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                              "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}

# Research-named outputs -> Ladrillo data/ file and the label transform tools/import_l27_data.py applies.
SHIPPED = [
    ("observations/calibration_targets.csv", "outputs/recalib_targets_ext.csv", "copy"),
    ("observations/calibration_targets_gsic_adj.csv", "outputs/recalib_targets_ext_gsicadj.csv", "copy"),
    ("drivers/glacier_block_temperature.csv", "data/observations/t_glac_blocks.csv", "header"),
    ("constants/glacier_block_constants.csv", "outputs/extc_block_constants.csv", "block_label"),
    ("drivers/greenland_zone_temperature.csv", "data/observations/t_gis_zones.csv", "copy"),
    ("calibration/gis_amp_prior.csv", "outputs/gis_amp_prior.csv", "copy"),
    ("constants/gis_amp_shape.csv", "outputs/gis_amp_shape.csv", "copy"),
    ("constants/gis_amp_shape_meta.csv", "outputs/gis_amp_shape_meta.csv", "copy"),
]
# Intermediates: research files Ladrillo does not ship; their research sha256 is in reference.sha256.
INTERMEDIATES = ["outputs/diag_constraint_anatomy_regions.csv", "outputs/lws_grace_extension_L14.csv",
                 "outputs/gis_driver_constants.csv", "outputs/gis_amp_anchor.csv"]
LADRILLO_OUT = "ladrillo"           # <work-dir>/ladrillo/<data-relative path>

# The chain. {S} = sources dir, {W} = work dir, {C} = Ladrillo commit. Order matters.
STEPS = [
    ("regions", ["build_region_inventory.py", "--sources", "{S}", "--work", "{W}"]),
    ("lws", ["build_lws_grace_extension.py", "--sources", "{S}", "--work", "{W}"]),
    ("targets", ["prep_recalib_targets_ext.py", "--sources", "{S}", "--work", "{W}", "--commit", "{C}"]),
    ("glacier", ["ladrillo_data.py", "--sources", "{S}", "--work", "{W}"]),
    ("greenland", ["build_t_gis.py", "--sources", "{S}", "--work", "{W}"]),
    ("amp_anchor", ["diag_gis_amp_anchor.py", "--sources", "{S}", "--work", "{W}"]),
    ("amp_shape", ["diag_gis_amp_cmip6.py", "--sources", "{S}", "--work", "{W}", "--commit", "{C}"]),
]
# Maintainer extractions (--from-raw). {R} = raw dir, {E} = extraction output dir (a copy of data/sources).
EXTRACT_STEPS = [
    ("x_hadcrut5", ["extract/build_t_glac.py", "--raw-dir", "{R}", "--sources", "{S0}", "--out-dir", "{E}/hadcrut5"]),
    ("x_greenland", ["extract/extract_greenland_zones.py", "--raw-dir", "{R}", "--out-dir", "{E}/greenland_zones"]),
    ("x_grace", ["extract/extract_grace_lws_domain.py", "--raw-dir", "{R}", "--sources", "{S0}",
                 "--out-dir", "{E}/grace"]),
    ("x_dangendorf", ["extract/extract_dangendorf_gmsl.py", "--raw-dir", "{R}", "--out-dir", "{E}/dangendorf2024"]),
]
EXTRACTED_FILES = ["hadcrut5/t_glac_hadcrut5.csv", "hadcrut5/t_glac_regions_hadcrut5.csv",
                   "greenland_zones/greenland_zone_series_HadCRUT5.csv",
                   "greenland_zones/greenland_zone_series_BerkeleyEarth.csv",
                   "greenland_zones/greenland_zone_series_GISTEMP.csv",
                   "grace/grace_lws_domain_monthly.csv", "dangendorf2024/dangendorf2024_gmsl_annual.csv"]
SEEDS = [dict(seed=2026, name="FIT_RNG_SEED", script="ladrillo_data.py",
              applied="four_rung_fit(): multi-start Nelder-Mead jitter; one stream per block via "
                      "block_rng(block) = default_rng([2026, sha256(name|basis|amp_b)[:8]])",
              affects="constants/glacier_block_constants.csv (b_fit_*, T_off_fit_*, kappa_anch_*, nu_anch_*)")]


def git(*args):
    r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True)
    return r.stdout.strip()


def check_inputs(root: Path, pins: dict, prefix: str, skip=()) -> list[str]:
    bad = []
    for rel, sha in sorted(pins.items()):
        if not rel.startswith(prefix) or rel[len(prefix):] in skip:
            continue
        p = root / rel[len(prefix):]
        if not p.exists():
            bad.append(f"MISSING  {rel}  (expected at {p})")
        elif TI.sha256(p) != sha:
            bad.append(f"MISMATCH {rel}  ({p})")
    return bad


def run(name, argv, snap: Path, subs: dict, log: list):
    cmd = [sys.executable, str(snap / argv[0])] + [a.format(**subs) for a in argv[1:]]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=str(snap), env={**os.environ, **BLAS_ENV}, capture_output=True, text=True)
    wall = time.time() - t0
    log.append(dict(step=name, cmd=" ".join(cmd), wall_s=round(wall, 2), returncode=r.returncode))
    tail = (r.stdout + r.stderr).strip().splitlines()[-3:]
    print(f"  [{name:12s}] {wall:6.1f} s  " + (" | ".join(tail) if r.returncode == 0 else "FAILED"))
    if r.returncode != 0:
        print(r.stdout[-4000:] + r.stderr[-4000:])
        sys.exit(f"step {name} failed (exit {r.returncode}); nothing installed")


def relabel(work: Path, snap: Path):
    sys.path.insert(0, str(snap))
    import import_l27_data as IMP   # the label-only transform that made data/ (SLOWP->SLOWG, FAST->FASTG)
    for dest, research_rel, how in SHIPPED:
        d = work / LADRILLO_OUT / dest
        d.parent.mkdir(parents=True, exist_ok=True)
        IMP.transform(str(work / research_rel), str(d), how)
        assert IMP.verify(str(work / research_rel), str(d), how), f"label transform not label-only: {dest}"


def versions():
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for m in ("numpy", "pandas", "scipy", "xarray", "netCDF4", "openpyxl", "shapefile", "matplotlib"):
        try:
            mod = __import__(m)
            out[m] = getattr(mod, "__version__", "?")
        except ImportError:
            out[m] = None
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--inputs-dir", required=True, help="the fetch-yourself inputs (fetch_inputs.py layout)")
    ap.add_argument("--sources", default=str(TI.SOURCES_DIR), help="the shipped inputs (default: data/sources/)")
    ap.add_argument("--from-raw", action="store_true",
                    help="re-run every maintainer extraction from <inputs-dir>/raw/ first")
    ap.add_argument("--tools-dir", default=str(HERE), help="scripts to snapshot (default: this directory)")
    ap.add_argument("--install", action="store_true", help="copy the eight built files into data/ (see docstring)")
    a = ap.parse_args()

    work = Path(a.work_dir).resolve()
    if work.exists() and any(work.iterdir()):
        sys.exit(f"--work-dir {work} is not empty; use a fresh directory")
    work.mkdir(parents=True, exist_ok=True)
    sources0 = Path(a.sources).resolve()
    inputs = Path(a.inputs_dir).expanduser().resolve()
    pins = TI.read_pins()
    commit = git("rev-parse", "--short", "HEAD") or "unknown"
    dirty = bool(git("status", "--porcelain", "--", "tools/targets", "data/sources"))
    print(f"build_targets | Ladrillo {commit}{' (+ uncommitted changes in tools/targets or data/sources)' if dirty else ''}"
          f" | work {work}")

    snap = work / "tools_snapshot"
    shutil.copytree(Path(a.tools_dir), snap, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copy2(REPO / "tools" / "import_l27_data.py", snap / "import_l27_data.py")
    log: list = []

    # every input checked BEFORE anything runs; extraction outputs are skipped when they will be re-made
    skip = set(EXTRACTED_FILES) if a.from_raw else set()
    bad = check_inputs(sources0, pins, "sources/", skip) + check_inputs(inputs, pins, "fetch/", skip)
    if a.from_raw:
        bad += check_inputs(inputs / "raw", pins, "raw/")
    if bad:
        print("\n".join("  " + b for b in bad))
        sys.exit("ABORT: inputs do not match inputs.sha256 (fetch_inputs.py --verify says more); nothing was run")
    n_in = sum(k.startswith(("sources/", "fetch/")) or (a.from_raw and k.startswith("raw/")) for k in pins) - len(skip)
    print(f"inputs: all {n_in} files match inputs.sha256")
    sources = work / "sources"                       # shipped + fetched, one tree
    shutil.copytree(sources0, sources)
    for k in pins:
        if k.startswith("fetch/") and k[len("fetch/"):] not in skip:
            rel = k[len("fetch/"):]
            (sources / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(inputs / rel, sources / rel)
    if a.from_raw:
        for rel in EXTRACTED_FILES:
            if (sources / rel).exists():
                (sources / rel).unlink()
        print(f"maintainer extractions from {inputs / 'raw'}:")
        for name, argv in EXTRACT_STEPS:
            run(name, argv, snap, dict(R=str(inputs / "raw"), S0=str(sources0), E=str(sources)), log)

    print("chain:")
    for name, argv in STEPS:
        run(name, argv, snap, dict(S=str(sources), W=str(work), C=commit), log)
    relabel(work, snap)

    outputs = {f"{LADRILLO_OUT}/{d}": TI.sha256(work / LADRILLO_OUT / d) for d, _, _ in SHIPPED}
    outputs.update({r: TI.sha256(work / r) for _, r, _ in SHIPPED})
    outputs.update({r: TI.sha256(work / r) for r in INTERMEDIATES})
    if a.from_raw:
        outputs.update({f"sources/{r}": TI.sha256(sources / r) for r in EXTRACTED_FILES})
    prov = dict(
        built=dt.datetime.now().isoformat(timespec="seconds"),
        ladrillo_commit=commit, uncommitted_changes=dirty,
        sources=str(sources0), inputs_dir=str(inputs), from_raw=a.from_raw,
        steps=log,
        seeds=SEEDS,
        randomness_note="No other step draws random numbers. The research diag_gis_amp_cmip6.py bootstraps "
                        "(seed 2026) never reached the shape and are not ported.",
        scripts={p.relative_to(snap).as_posix(): TI.sha256(p) for p in sorted(snap.rglob("*.py"))},
        inputs={k: v for k, v in pins.items() if k.startswith(("sources/", "fetch/")) or
                (a.from_raw and k.startswith("raw/"))},
        outputs=outputs,
        versions=versions(),
        blas_env=BLAS_ENV,
    )
    (work / "BUILD_PROVENANCE.json").write_text(json.dumps(prov, indent=1) + "\n")
    print(f"built {len(SHIPPED)} Ladrillo files under {work / LADRILLO_OUT} (+ {len(INTERMEDIATES)} intermediates); "
          f"wrote BUILD_PROVENANCE.json")

    if a.install:
        for d, _, _ in SHIPPED:
            shutil.copy2(work / LADRILLO_OUT / d, REPO / "data" / d)
        print("INSTALLED into data/. data/MANIFEST.csv's sha256 column now disagrees for any file that changed; "
              "check_targets_identity.py will say which.")


if __name__ == "__main__":
    main()
