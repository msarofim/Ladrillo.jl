#!/usr/bin/env python3
"""
mutation_test_tables.py -- mutation test of check_tables_identity.py: break what the gate guards, on a SCRATCH COPY of
python/tables/ and data/{constants,posterior}/, and confirm the gate FAILS.

  python python/tables/mutation_test_tables.py <SLR-RFF-BRICK repo> <scratch dir>

Each mutant is a fresh copy with ONE change; the four table scripts run on it and the gate compares its outputs with
the research outputs. Nothing in this repository is modified. Table 4 reads the hindcast bands in outputs/hindcast/
(scripts/posterior_predictive.jl and scripts/posterior_predictive_brick20.jl), which must exist; they are copied too. The run that certified the port is saved in
python/tables/mutation_result.txt.

THE WITNESS. Whether a mutant's outputs changed at all is decided independently of the gate: its five output files
are compared BYTE FOR BYTE with those of an unmutated copy run the same way, after blanking dates and times (the only
run-dependent text; both copies sit outside git, so both provenance commits are empty). The gate must FAIL every
mutant whose outputs changed and PASS every mutant whose outputs did not (and the unmutated copy itself).

  mutant              change                                                       expected
  posterior_1ulp      one posterior value (ais_slope, draw 1) moved by 1 ulp       outputs unchanged, gate PASS:
                                                                                   pandas' default float parser (which
                                                                                   the research scripts use) reads
                                                                                   both texts as the same double
  posterior_1e-12     the same value moved by a relative 1e-12                     FAIL
  geoprior_sd_1ulp    the paleo-prior sd of ais_slope moved by 1 ulp               FAIL
  prop_min_0.78       Table A2's PROP_MIN 0.80 -> 0.78 (PC 3's ratio 0.789 crosses) FAIL
  prop_min_0.85       Table A2's PROP_MIN 0.80 -> 0.85 (no ratio crosses; the       FAIL
                      caption states the threshold)
  double_rename       gic_a_FASTG -> gic_a_FASTGG in the priors, the posterior      FAIL
                      header and the Table A1 script (a rename applied twice)
  ln_prior_unscaled   Table A1 md: the ln display scales the posterior quantiles     FAIL
                      but not the prior (the prior column stays log10)
  ln10_2.30           Table A1: ln 10 written as 2.30                                SCRIPT: [LN-DISPLAY] stops it
  t4_window           Table 4: the 1950-1992 window ends in 1993                     FAIL
  t4_band             Table 4: Ladrillo's band read from the predictive p05          FAIL
  t4_obs_1e-12        Table 4: one Ladrillo-file glacier observation moved 1e-12     FAIL
  t4_brick_1e-12      Table 4: one BRICK 2.0 Greenland median moved 1e-12            FAIL
  t4_refwin           Table 4: re-reference window 1995-2006 (targets are zeroed on   SCRIPT: the zeroing check
                      1995-2005)                                                     stops it
SCRIPT means the mutated script must refuse to write (non-zero exit), so there is no output to gate.
"""
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
RES = os.path.abspath(sys.argv[1])
SCRATCH = os.path.abspath(sys.argv[2])
SCRIPTS = ("diag_ais_block_pca.py", "ladrillo_table_a2.py", "ladrillo_prior_posterior_table.py",
           "ladrillo_vs_brick20_scorecard.py")
OUTPUTS = ("diag_ais_block_pca.csv", "ladrillo_table_a2.md", "ladrillo_prior_posterior.csv", "ladrillo_prior_posterior.md",
           "ladrillo_vs_brick20_scorecard.csv")
POST = os.path.join("data", "posterior", "ladrillo_L27_posterior.csv")
PRIORS = os.path.join("data", "constants", "ladrillo_priors_L27.csv")
GEO = os.path.join("data", "constants", "paleo_geo_prior_ton.csv")
A2 = os.path.join("python", "tables", "ladrillo_table_a2.py")
A1 = os.path.join("python", "tables", "ladrillo_prior_posterior_table.py")
T4 = os.path.join("python", "tables", "ladrillo_vs_brick20_scorecard.py")
HIND_L = os.path.join("outputs", "hindcast", "postpred_components_timeseries.csv")
HIND_B = os.path.join("outputs", "hindcast", "postpred_brick20_components_timeseries.csv")
TARGETS = os.path.join("data", "observations", "calibration_targets.csv")
THRESHOLDS = "IDENT_MAX, PROP_MIN, LOAD_MIN = 0.20, 0.80, 0.30"
DATETIME_RX = re.compile(r"\d{4}-\d{2}-\d{2}( \d{2}:\d{2})?")


def edit(root, rel, fn):
    p = os.path.join(root, rel)
    s = open(p, encoding="utf-8").read()
    t = fn(s)
    assert t != s, f"mutation of {rel} changed nothing"
    open(p, "w", encoding="utf-8").write(t)


def move_cell(s, line_no, col_name, how):
    """Move one CSV cell (data line line_no, 1-based; column col_name) by `how`: 'ulp' or a relative step."""
    lines = s.split("\n")
    j = lines[0].split(",").index(col_name)
    cells = lines[line_no].split(",")
    x = float(cells[j])
    cells[j] = repr(float(np.nextafter(x, np.inf)) if how == "ulp" else x * (1 + how))
    lines[line_no] = ",".join(cells)
    return "\n".join(lines)


def ulp_geo_sd(s, name):
    lines = s.split("\n")
    j = next(ln for ln in lines if ln.startswith("names,")).split(",").index(name)
    k = next(i for i, ln in enumerate(lines) if ln.startswith("sd,"))
    cells = lines[k].split(",")
    cells[j] = repr(float(np.nextafter(float(cells[j]), np.inf)))
    lines[k] = ",".join(cells)
    return "\n".join(lines)


MUTANTS = [
    ("posterior_1ulp", "PASS", lambda r: edit(r, POST, lambda s: move_cell(s, 1, "ais_slope", "ulp"))),
    ("posterior_1e-12", "FAIL", lambda r: edit(r, POST, lambda s: move_cell(s, 1, "ais_slope", 1e-12))),
    ("geoprior_sd_1ulp", "FAIL", lambda r: edit(r, GEO, lambda s: ulp_geo_sd(s, "ais_slope"))),
    ("prop_min_0.78", "FAIL", lambda r: edit(r, A2, lambda s: s.replace(THRESHOLDS, THRESHOLDS.replace("0.80", "0.78")))),
    ("prop_min_0.85", "FAIL", lambda r: edit(r, A2, lambda s: s.replace(THRESHOLDS, THRESHOLDS.replace("0.80", "0.85")))),
    ("double_rename", "FAIL", lambda r: (edit(r, PRIORS, lambda s: s.replace("gic_a_FASTG,", "gic_a_FASTGG,")),
                                         edit(r, POST, lambda s: s.replace("gic_a_FASTG,", "gic_a_FASTGG,", 1)),
                                         edit(r, A1, lambda s: s.replace('"gic_a_FASTG"', '"gic_a_FASTGG"')))),
    ("ln_prior_unscaled", "FAIL", lambda r: edit(r, A1, lambda s: s.replace("                    pr[k] = pr[k] * LN10\n",
                                                                            "                    pass\n"))),
    ("ln10_2.30", "SCRIPT", lambda r: edit(r, A1, lambda s: s.replace("LN10 = np.log(10.0)", "LN10 = 2.30"))),
    ("t4_window", "FAIL", lambda r: edit(r, T4, lambda s: s.replace("(1950, 1992)", "(1950, 1993)"))),
    ("t4_band", "FAIL", lambda r: edit(r, T4, lambda s: s.replace('{"Ladrillo": "p05",', '{"Ladrillo": "pred_p05",'))),
    ("t4_obs_1e-12", "FAIL", lambda r: edit(r, HIND_L, lambda s: move_cell(s, 60, "glaciers_obs", 1e-12))),
    ("t4_brick_1e-12", "FAIL", lambda r: edit(r, HIND_B, lambda s: move_cell(s, 110, "gis_p50", 1e-12))),
    ("t4_refwin", "SCRIPT", lambda r: edit(r, T4, lambda s: s.replace("REF = (1995, 2005)", "REF = (1995, 2006)"))),
]


def fresh(name):
    root = os.path.join(SCRATCH, name)
    shutil.rmtree(root, ignore_errors=True)
    for d in (os.path.join("python", "tables"), os.path.join("data", "constants"), os.path.join("data", "posterior")):
        shutil.copytree(os.path.join(REPO, d), os.path.join(root, d), ignore=shutil.ignore_patterns("__pycache__"))
    for f in (HIND_L, HIND_B, TARGETS):
        if not os.path.isfile(os.path.join(REPO, f)):
            sys.exit(f"{f} is missing: run scripts/posterior_predictive.jl and scripts/posterior_predictive_brick20.jl")
        os.makedirs(os.path.dirname(os.path.join(root, f)), exist_ok=True)
        shutil.copy2(os.path.join(REPO, f), os.path.join(root, f))
    return root


def outputs(root):
    return {f: DATETIME_RX.sub("<date>", open(os.path.join(root, "outputs", "tables", f), encoding="utf-8").read())
            for f in OUTPUTS}


def run(root):
    out = os.path.join(root, "outputs", "tables")
    for s in SCRIPTS:
        p = subprocess.run([sys.executable, os.path.join(root, "python", "tables", s), f"--outdir={out}"],
                           capture_output=True, text=True)
        if p.returncode:
            return "SCRIPT", 0, [f"{s} stopped: {p.stderr.strip().splitlines()[-1]}"]
    g = subprocess.run([sys.executable, os.path.join(root, "python", "tables", "check_tables_identity.py"), RES,
                        f"--pkgdir={out}"], capture_output=True, text=True)
    summary = [ln.strip() for ln in g.stdout.splitlines()
               if ln.strip().startswith(("IDENTICAL", "FAIL", "MISSING")) or "bit-identical" in ln]
    ndiff = sum(ln.strip().startswith("DIFF") for ln in g.stdout.splitlines())
    return ("PASS" if g.returncode == 0 else "FAIL"), ndiff, summary


ctrl = fresh("control")
cv, cn, csumm = run(ctrl)
base = outputs(ctrl)
rows = [("control (unmutated copy)", "PASS", cv, cn, [], csumm)]
bad = cv != "PASS"
for name, expect, mutate in MUTANTS:
    root = fresh(name)
    mutate(root)
    v, n, summ = run(root)
    if v == "SCRIPT":       ## the mutated script refused to write: nothing to gate, and nothing changed on disk
        changed, consistent = [], v == expect
    else:
        changed = [f for f, t in outputs(root).items() if t != base[f]]
        ## the gate must fail exactly when an output changed; and the expected change must have happened
        consistent = (v == "FAIL") == bool(changed) and v == expect
    bad += not consistent
    rows.append((name, expect, v, n, changed, summ))
for name, expect, v, n, changed, summ in rows:
    print(f"{name:26s} expected {expect}  gate {v} ({n:3d} DIFF lines)  outputs changed: "
          f"{', '.join(changed) if changed else 'none'}  "
          f"{'OK' if v == expect and (v != 'PASS' or not changed) and (v != 'FAIL' or changed) else 'UNEXPECTED'}")
    for s in summ:
        print(f"{'':30s}{s}")
print(f"[TABLES-MUTATION] {'PASS: every mutant behaved as expected' if not bad else f'FAIL: {bad} unexpected'}")
sys.exit(1 if bad else 0)
