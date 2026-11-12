#!/usr/bin/env python3
"""
mutation_test_targets.py -- does the targets identity gate catch what it guards?

Each mutant copies tools/targets/ (or an input tree), changes ONE thing,
rebuilds with build_targets.py --tools-dir <copy>, and runs
check_targets_identity.py on the result. A builder mutant must make the gate
FAIL and name the files it should; an input mutant must make build_targets.py
abort before any step runs; the controls must PASS (a change that cannot reach
any output) or fail ONLY on the intermediate it touches (showing what the
intermediate rows add).

    python tools/targets/mutation_test_targets.py --inputs-dir ~/ladrillo_inputs [--out result.txt]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
import target_inputs as TI  # noqa: E402

# (name, file under tools/targets, old text, new text, expectation, files that must FAIL (substring match))
BUILDER_MUTANTS = [
    ("M1 altimetry sigma 4.0 -> 4.1 mm", "prep_recalib_targets_ext.py",
     "ALT_SIGMA_MM     = 4.0 ", "ALT_SIGMA_MM     = 4.1 ", "FAIL",
     ["calibration_targets.csv", "recalib_targets_ext.csv"]),
    ("M2 land-water read-back emulation removed (1-ulp)", "build_lws_grace_extension.py",
     "fred_lws = as_read_back_from_targets_csv(frederikse_tws(a.sources))", "fred_lws = frederikse_tws(a.sources)",
     "FAIL", ["lws_grace_extension_L14.csv", "calibration_targets.csv"]),
    ("M3 glacier fit seed 2026 -> 2027", "ladrillo_data.py",
     "FIT_RNG_SEED = 2026\n", "FIT_RNG_SEED = 2027\n", "FAIL", ["glacier_block_constants.csv"]),
    ("M4 amp shape support to 3.25 K (the fullcurve arm)", "diag_gis_amp_cmip6.py",
     "SHAPE_DT_MAX = 2.75", "SHAPE_DT_MAX = 3.25", "FAIL", ["gis_amp_shape.csv", "gis_amp_shape_meta.csv"]),
    ("M5 Greenland driver from Berkeley Earth instead of HadCRUT5", "build_t_gis.py",
     'DRIVER_PRODUCT = "HadCRUT5"', 'DRIVER_PRODUCT = "BerkeleyEarth"', "FAIL", ["greenland_zone_temperature.csv"]),
    ("C1 control: a console message only", "prep_recalib_targets_ext.py",
     'print(f"\\nWrote {OUT_CSV}', 'print(f"\\nWROTE {OUT_CSV}', "PASS", []),
    ("C2 control: melt-rate smoothing 11 -> 9 yr (feeds only an intermediate)", "build_t_gis.py",
     "SMOOTH = 11 ", "SMOOTH = 9 ", "FAIL-ONLY", ["gis_driver_constants.csv"]),
]
INPUT_MUTANTS = [
    ("I1 one digit in a shipped input (GlaMBIE global)", "sources", "glambie/calendar_years/0_global.csv"),
    ("I2 one byte in a fetched input (Frederikse workbook)", "fetch",
     "frederikse2020/frederikse2020_global_basin_timeseries.xlsx"),
    ("I3 one digit in a shipped rolling-product file (NOAA STAR annual means)", "sources", "noaa/nasa_gmsl_annual.csv"),
]


def flip_digit(p: Path):
    b = bytearray(p.read_bytes())
    i = b.rfind(b"5")
    assert i > 0, f"no digit 5 in {p}"
    b[i] = ord("6")
    p.write_bytes(bytes(b))


def build(work, tools, inputs, sources=None):
    cmd = [sys.executable, str(HERE / "build_targets.py"), "--work-dir", str(work), "--inputs-dir", str(inputs),
           "--tools-dir", str(tools)] + (["--sources", str(sources)] if sources else [])
    return subprocess.run(cmd, capture_output=True, text=True)


def gate(work):
    r = subprocess.run([sys.executable, str(HERE / "check_targets_identity.py"), "--work-dir", str(work)],
                       capture_output=True, text=True)
    fails = [ln.split()[1] for ln in r.stdout.splitlines() if " FAIL " in f" {ln} " and ln.startswith("  ")]
    return r.returncode, fails, r.stdout.strip().splitlines()[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs-dir", required=True)
    ap.add_argument("--out", default=str(HERE / "mutation_result.txt"))
    a = ap.parse_args()
    inputs = Path(a.inputs_dir).expanduser().resolve()
    tmp = Path(tempfile.mkdtemp(prefix="targets_mutation_"))
    lines, ok_all = [], True
    commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], capture_output=True,
                            text=True).stdout.strip()
    lines.append(f"targets gate mutation test, {time.strftime('%Y-%m-%d %H:%M')}, Ladrillo {commit} "
                 f"(+ working tree), inputs {inputs}")
    try:
        for i, (name, fname, old, new, expect, must) in enumerate(BUILDER_MUTANTS):
            tools = tmp / f"tools_{i}"
            shutil.copytree(HERE, tools, ignore=shutil.ignore_patterns("__pycache__"))
            src = (tools / fname).read_text()
            assert src.count(old) == 1, f"{name}: mutation site found {src.count(old)} times"
            (tools / fname).write_text(src.replace(old, new))
            work = tmp / f"work_{i}"
            b = build(work, tools, inputs)
            if b.returncode != 0:
                verdict, detail = "BUILD-ERROR", (b.stdout + b.stderr).strip().splitlines()[-1]
                ok = False
            else:
                rc, fails, summary = gate(work)
                if expect == "PASS":
                    ok = rc == 0
                elif expect == "FAIL":
                    ok = rc != 0 and all(any(m in f for f in fails) for m in must)
                else:   # FAIL-ONLY: fails, and only on the named files
                    ok = rc != 0 and fails and all(any(m in f for m in must) for f in fails)
                verdict = "PASS" if rc == 0 else "FAIL"
                detail = f"{summary}; failing: {', '.join(fails) if fails else '-'}"
            ok_all &= ok
            lines.append(f"{'ok ' if ok else 'BAD'} {name:62s} expect {expect:9s} got {verdict:11s} {detail}")
            print(lines[-1], flush=True)
        for j, (name, kind, rel) in enumerate(INPUT_MUTANTS):
            src_tree = tmp / f"sources_{j}"
            in_tree = tmp / f"inputs_{j}"
            shutil.copytree(TI.SOURCES_DIR, src_tree)
            shutil.copytree(inputs, in_tree, symlinks=True)
            flip_digit((src_tree if kind == "sources" else in_tree) / rel)
            work = tmp / f"work_in_{j}"
            b = build(work, HERE, in_tree, sources=src_tree)
            out = b.stdout + b.stderr
            ok = b.returncode != 0 and "nothing was run" in out and "chain:" not in out and f"MISMATCH {kind}/{rel}" in out
            ok_all &= ok
            lines.append(f"{'ok ' if ok else 'BAD'} {name:62s} expect ABORT     got "
                         f"{'ABORT' if b.returncode else 'RAN':11s} {out.strip().splitlines()[-1]}")
            print(lines[-1], flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    lines.append("MUTATION TEST " + ("PASS: every mutant behaved as expected" if ok_all else "FAIL"))
    print(lines[-1])
    Path(a.out).write_text("\n".join(lines) + "\n")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
