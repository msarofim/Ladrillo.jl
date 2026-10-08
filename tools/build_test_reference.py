#!/usr/bin/env python3
"""
build_test_reference.py -- extract the test suite's reference values from the RESEARCH code's outputs.

    python tools/build_test_reference.py <SLR-RFF-BRICK repo>

Maintainer-only: run once per certified release. test/runtests.jl then needs only test/reference/, so a reviewer runs the
suite from a fresh clone. The references come from the research drivers (an independent implementation of the same
model), not from this package, so the suite keeps testing this package against the code that produced the paper.

Lines are COPIED as text, never parsed and re-printed, so every value is exactly what the research driver wrote.
Writes test/reference/*.csv and test/reference/SOURCES.csv (file, source, research commit, source sha256).
"""
import hashlib, os, subprocess, sys

RES = os.path.abspath(sys.argv[1])
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "test", "reference")
TAP = "_tap4p69K_V5p64m_tau800"
## Joint-arm draws to re-run (index into posterior rows 1:5:10000). Evenly spaced, plus the three draws whose fast
## dynamics switches on inside the record on their own FaIR config (CHANGELOG 2026-10-08c), so the threshold path is
## exercised too.
REGRESSION_DRAWS = sorted(set(range(1, 2001, 67)) | {229, 610, 1260})
SCENARIOS = ("ssp245", "ssp585")
## (reference file, research source, kind): kind "draws" keeps only REGRESSION_DRAWS, "copy" copies the whole file
FILES = [(f"research_joint_draws_{s}.csv", f"outputs/scope_slr_fairunc_draws_{s}_spliced_L27{TAP}.csv", "draws")
         for s in SCENARIOS] + \
        [("research_brick20_joint_draws_ssp245.csv", "outputs/scope_slr_fairunc_draws_ssp245_spliced_oldbrick.csv",
          "draws")] + \
        [(f"research_postpred_{stem}.csv", f"outputs/postpred_L27_{stem}.csv", "copy")
         for stem in ("components_timeseries", "bias", "coverage")]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


os.makedirs(OUT, exist_ok=True)
commit = subprocess.run(["git", "-C", RES, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
dirty = subprocess.run(["git", "-C", RES, "status", "--porcelain", "--"] + [s for _, s, _ in FILES],
                       capture_output=True, text=True).stdout.strip()
if dirty:
    sys.exit(f"research sources have uncommitted changes, refusing:\n{dirty}")
rows = ["file,source,research_commit,source_sha256,lines_kept,regression_draws"]
for name, src, kind in FILES:
    sp = os.path.join(RES, src)
    with open(sp) as f:
        header = f.readline()
        if kind == "draws":
            cols = header.rstrip("\n").split(",")
            assert cols[0] == "draw", (src, cols[:3])
            keep = [ln for ln in f if int(ln.split(",", 1)[0]) in REGRESSION_DRAWS]
        else:
            keep = f.readlines()
    with open(os.path.join(OUT, name), "w") as g:
        g.write(header); g.writelines(keep)
    rows.append(f"{name},{src},{commit},{sha256(sp)},{len(keep)},"
                f"{' '.join(map(str, REGRESSION_DRAWS)) if kind == 'draws' else ''}")
    print(f"  {name}: {len(keep)} lines from {src}")
with open(os.path.join(OUT, "SOURCES.csv"), "w") as g:
    g.write("\n".join(rows) + "\n")
print(f"wrote {os.path.relpath(OUT, HERE)} from SLR-RFF-BRICK {commit}")
