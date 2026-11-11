#!/usr/bin/env python3
"""
build_test_reference.py -- extract the test suite's reference values from the RESEARCH code's outputs.

    python tools/build_test_reference.py <SLR-RFF-BRICK repo>
    python tools/build_test_reference.py <SLR-RFF-BRICK repo> --only research_ladrillo_table_a2_L27_sub10k.md ...

Maintainer-only: run once per certified release. test/runtests.jl then needs only test/reference/, so a reviewer runs the
suite from a fresh clone. The references come from the research drivers (an independent implementation of the same
model), not from this package, so the suite keeps testing this package against the code that produced the paper.

Lines are COPIED as text, never parsed and re-printed, so every value is exactly what the research driver wrote.
Writes test/reference/*.csv (and the Table A2 *.md) and test/reference/SOURCES.csv (file, source, research commit,
source sha256). --only copies just the named reference files and rewrites (or appends) only their SOURCES.csv rows; every
other row is kept byte-for-byte, because its file was copied at its own research commit.
"""
import argparse, hashlib, os, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("research", help="the SLR-RFF-BRICK repository")
ap.add_argument("--only", nargs="+", metavar="FILE", help="copy only these reference files (names under test/reference)")
ARGS = ap.parse_args()
RES = os.path.abspath(ARGS.research)
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
         for stem in ("components_timeseries", "bias", "coverage")] + \
        [("research_ladrillo_table_a2_L27_sub10k.md", "outputs/ladrillo_table_a2_L27_sub10k.md", "copy"),   # Table A2
         ("research_ladrillo_prior_posterior_L27.csv", "outputs/ladrillo_prior_posterior_L27.csv", "copy"),   # Table A1
         ("research_ladrillo_prior_posterior_L27.md", "outputs/ladrillo_prior_posterior_L27.md", "copy"),     # Table A1
         ("research_diag_ais_block_pca_L27_sub10k.csv", "outputs/diag_ais_block_pca_L27_sub10k.csv", "copy"),  # A2 basis
         ("research_postpred_brick20_components_timeseries.csv", "outputs/postpred_oldbrick_components_timeseries.csv",
          "copy"),                                                                                          # Table 4 input
         ("research_scope_ladrillo_vs_brick20_scorecard_L27.csv", "outputs/scope_ladrillo_vs_brick20_scorecard_L27.csv",
          "copy")]                                                                                          # Table 4


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


if ARGS.only:
    unknown = sorted(set(ARGS.only) - {n for n, _, _ in FILES})
    if unknown:
        sys.exit(f"--only: not in FILES: {unknown}")
TODO = [f for f in FILES if not ARGS.only or f[0] in ARGS.only]
os.makedirs(OUT, exist_ok=True)
commit = subprocess.run(["git", "-C", RES, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
dirty = subprocess.run(["git", "-C", RES, "status", "--porcelain", "--"] + [s for _, s, _ in TODO],
                       capture_output=True, text=True).stdout.strip()
if dirty:
    sys.exit(f"research sources have uncommitted changes, refusing:\n{dirty}")
HEADER = "file,source,research_commit,source_sha256,lines_kept,regression_draws"
rows = [HEADER]
for name, src, kind in TODO:
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
    if kind == "copy":      # a whole-file copy must be byte-identical (no newline or encoding translation)
        assert sha256(os.path.join(OUT, name)) == sha256(sp), f"{name}: copy differs from {src}"
    rows.append(f"{name},{src},{commit},{sha256(sp)},{len(keep)},"
                f"{' '.join(map(str, REGRESSION_DRAWS)) if kind == 'draws' else ''}")
    print(f"  {name}: {len(keep)} lines from {src}")
if ARGS.only:
    ## keep every other row as it is; replace a copied file's row in place, append a new one at the end
    old = open(os.path.join(OUT, "SOURCES.csv")).read().split("\n")
    assert old[0] == HEADER and old[-1] == "", "SOURCES.csv: unexpected header or no final line break"
    old = old[:-1]
    at = {ln.split(",", 1)[0]: i for i, ln in enumerate(old) if i > 0}
    for r in rows[1:]:
        if r.split(",", 1)[0] in at:
            old[at[r.split(",", 1)[0]]] = r
        else:
            old.append(r)
    rows = old
with open(os.path.join(OUT, "SOURCES.csv"), "w") as g:
    g.write("\n".join(rows) + "\n")
print(f"wrote {os.path.relpath(OUT, HERE)} from SLR-RFF-BRICK {commit}")
