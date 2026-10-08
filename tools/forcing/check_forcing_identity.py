"""
check_forcing_identity.py -- the IDENTITY GATE: does a rebuild reproduce every
shipped forcing file EXACTLY?

For each file Ladrillo ships in data/forcing/ (enumerated from data/MANIFEST.csv,
whose recorded sha256 is checked against the reference copy first, so the gate
cannot be passed by a reference that was edited in place):

  1. sha256 of the candidate == sha256 of data/forcing/<file>     -> PASS (bytes)
  2. otherwise, read both with pandas `float_precision="round_trip"` (so the
     parse is exact: every printed float maps to the one float64 that printed
     it) and require identical column names, identical shape, and EVERY CELL
     `==` (NaN equal only to NaN in the same cell; text columns compared as
     text)                                                         -> PASS (values)
  3. anything else                                                 -> FAIL, with the
     max |difference| over numeric cells and where it is.

There is NO tolerance. A file that differs in its last bit FAILS and the report
says by how much; that number is the diagnosis, not a reason to pass it.

    python check_forcing_identity.py --candidate-dir WORK/forcing
    python check_forcing_identity.py --candidate-dir WORK/forcing --files fair_mean_gmst_ssp245harm.csv
    python check_forcing_identity.py --candidate-dir WORK/forcing --allow-missing
    python check_forcing_identity.py --self-test       # mutation-test the gate itself

--allow-missing gates a PARTIAL build (build_forcing.py skips what needs an
absent input): every candidate file that EXISTS must pass exactly as above; the
expected files that do not exist are listed as ABSENT and do not fail the gate.
At least one file must exist. Without the flag a missing file is a FAIL.

Secondary mode, for diagnosis in the maintainer's research checkout: compare the
INTERMEDIATES a build wrote (calibration dirs, emissions files, van Vuuren
splices) with the research repo's copies:
    python check_forcing_identity.py --intermediates WORK --fairtofredi-root <research-checkout>
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from forcing_inputs import sha256  # noqa: E402

# -- Named constants ----------------------------------------------------------
LADRILLO_ROOT = HERE.parents[1]
REFERENCE_DIR = LADRILLO_ROOT / "data" / "forcing"
MANIFEST = LADRILLO_ROOT / "data" / "MANIFEST.csv"
MANIFEST_PREFIX = "forcing/"
READ_KW = dict(float_precision="round_trip")
SELF_TEST_FILE = "fair_mean_gmst_ssp245harm.csv"

# build_forcing.py work-dir subdir  <->  research-repo intermediate dir, files
# (None = every csv built), and the first year column compared (None = the whole
# file). The van Vuuren splice is built from vv_template_skeleton.csv, so its
# history columns hold Zenodo's own history where the research copy holds Smith
# 2024; nothing downstream reads them, so only the 2023.5+ columns (and the
# names) are compared there.
VV_SCENARIO_FIRST_COL = 2023.5
INTERMEDIATE_MAP = [
    ("calibration_v160", "calibration_v160",
     ["calibrated_constrained_parameters_1.6.0.csv",
      "species_configs_properties_1.6.0.csv", "volcanic_solar.csv",
      "emissions_v160_cmip7_ssp245_harmonized.csv"], None),
    ("calibration_v160_prod", "calibration_v160_prod", None, None),
    ("vanvuuren/spliced_ext_zenodo", "data/vanvuuren/spliced_ext_zenodo", None,
     VV_SCENARIO_FIRST_COL),
]


def compare(cand: Path, ref: Path):
    """-> dict(sha, values, max_abs, where, note). Never applies a tolerance."""
    r = {"sha": False, "values": False, "max_abs": None, "where": "", "note": ""}
    if not cand.exists():
        r["note"] = "candidate MISSING"
        return r
    if sha256(cand) == sha256(ref):
        r.update(sha=True, values=True, max_abs=0.0)
        return r
    return compare_frames(pd.read_csv(cand, **READ_KW), pd.read_csv(ref, **READ_KW), r)


def compare_frames(a, b, r=None):
    """Cell-by-cell comparison of two parsed tables (no tolerance)."""
    r = r or {"sha": False, "values": False, "max_abs": None, "where": "", "note": ""}
    if list(a.columns) != list(b.columns):
        r["note"] = (f"columns differ ({len(a.columns)} vs {len(b.columns)}; "
                     f"first mismatch "
                     f"{next((x, y) for x, y in zip(list(a.columns) + [None]*9, list(b.columns) + [None]*9) if x != y)})")
        return r
    if a.shape != b.shape:
        r["note"] = f"shape {a.shape} vs {b.shape}"
        return r
    all_eq, worst, where = True, 0.0, ""
    for c in a.columns:
        x, y = a[c], b[c]
        num = pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y)
        if num:
            xv, yv = x.to_numpy(float), y.to_numpy(float)
            same = (xv == yv) | (np.isnan(xv) & np.isnan(yv))
            if not same.all():
                all_eq = False
                d = np.abs(xv - yv)
                d[same] = 0.0
                i = int(np.nanargmax(d))
                if d[i] > worst or not where:
                    worst, where = float(d[i]), f"col {c!r} row {i}"
        else:
            if not x.astype(str).equals(y.astype(str)):
                all_eq = False
                where = where or f"text col {c!r}"
    r.update(values=all_eq, max_abs=worst, where=where)
    if all_eq:
        r["note"] = "bytes differ, every value identical (formatting only)"
    return r


def expected_files(manifest=MANIFEST):
    rows = list(csv.DictReader(open(manifest)))
    return {row["file"][len(MANIFEST_PREFIX):]: row["sha256"]
            for row in rows if row["file"].startswith(MANIFEST_PREFIX)}


def gate(cand_dir: Path, ref_dir: Path, files=None, report=None, allow_missing=False):
    exp = expected_files()
    names = files or sorted(exp)
    unknown = [n for n in names if n not in exp]
    if unknown:
        sys.exit(f"ERROR: not shipped forcing files (absent from MANIFEST): {unknown}")
    absent = [n for n in names if not (cand_dir / n).exists()] if allow_missing else []
    if allow_missing and len(absent) == len(names):
        sys.exit(f"ERROR: --allow-missing, but none of the {len(names)} expected "
                 f"files is in {cand_dir}; nothing to gate.")
    print(f"IDENTITY GATE: {cand_dir}  vs  {ref_dir}")
    print(f"  {len(names)} file(s); reference sha256 checked against "
          f"{MANIFEST.relative_to(LADRILLO_ROOT)} first"
          + (f"; --allow-missing: {len(absent)} ABSENT, not gated" if allow_missing else "")
          + "\n")
    hdr = f"{'file':<38} {'sha256':<7} {'values':<7} {'max|diff|':>11}  verdict"
    print(hdr)
    print("-" * len(hdr))
    rows, n_fail = [], 0
    for n in names:
        ref = ref_dir / n
        if sha256(ref) != exp[n]:
            sys.exit(f"ERROR: reference {ref} does not match MANIFEST.csv's sha256 "
                     f"-- the reference itself has changed; refusing to gate.")
        if n in absent:
            print(f"{n:<38} {'-':<7} {'-':<7} {'':>11}  ABSENT (not built)")
            rows.append({"file": n, "sha256_match": None, "value_match": None,
                         "max_abs_diff": None, "where": "", "note": "not built",
                         "verdict": "ABSENT"})
            continue
        r = compare(cand_dir / n, ref)
        ok = r["sha"] or r["values"]
        n_fail += not ok
        md = "" if r["max_abs"] is None else f"{r['max_abs']:.3e}"
        verdict = "PASS" if ok else "FAIL"
        extra = f"  [{r['note'] or r['where']}]" if (r["note"] or not ok) else ""
        print(f"{n:<38} {('yes' if r['sha'] else 'no'):<7} "
              f"{('yes' if r['values'] else 'no'):<7} {md:>11}  {verdict}{extra}")
        rows.append({"file": n, "sha256_match": r["sha"], "value_match": r["values"],
                     "max_abs_diff": r["max_abs"], "where": r["where"],
                     "note": r["note"], "verdict": verdict})
    print("-" * len(hdr))
    print(f"{len(names) - len(absent) - n_fail} PASS, {n_fail} FAIL"
          + (f", {len(absent)} ABSENT (not built; --allow-missing)" if absent else ""))
    if absent:
        print("ABSENT: " + " ".join(absent))
    if report:
        pd.DataFrame(rows).to_csv(report, index=False)
        print(f"report: {report}")
    return n_fail


def self_test():
    """Mutation-test the gate: it must FAIL a one-ulp change and a dropped row,
    and must PASS (values) a formatting-only change while FAILING sha256."""
    src = REFERENCE_DIR / SELF_TEST_FILE
    df = pd.read_csv(src, **READ_KW)
    col = df.columns[1]
    tmp = Path(tempfile.mkdtemp(prefix="forcing_gate_selftest_"))
    results = []
    try:
        # (a) identical copy -> PASS on sha
        shutil.copy2(src, tmp / "a.csv")
        r = compare(tmp / "a.csv", src)
        results.append(("identical copy", r["sha"] and r["values"], "PASS sha"))
        # (b) one cell moved by ONE ULP -> FAIL
        m = df.copy()
        i = len(m) // 2
        m.loc[i, col] = np.nextafter(m.loc[i, col], np.inf)
        m.to_csv(tmp / "b.csv", index=False)
        r = compare(tmp / "b.csv", src)
        results.append(("one cell +1 ulp", (not r["sha"]) and (not r["values"])
                        and r["max_abs"] > 0, f"FAIL, max|diff| {r['max_abs']:.3e}"))
        # (c) last row dropped -> FAIL
        df.iloc[:-1].to_csv(tmp / "c.csv", index=False)
        r = compare(tmp / "c.csv", src)
        results.append(("last row dropped", not (r["sha"] or r["values"]), r["note"]))
        # (d) formatting only (CRLF line endings) -> sha FAIL, values PASS
        (tmp / "d.csv").write_bytes(src.read_bytes().replace(b"\n", b"\r\n"))
        r = compare(tmp / "d.csv", src)
        results.append(("CRLF line endings", (not r["sha"]) and r["values"],
                        "sha no / values yes"))
    finally:
        shutil.rmtree(tmp)
    print(f"GATE SELF-TEST on {SELF_TEST_FILE}:")
    for name, ok, what in results:
        print(f"  {'ok ' if ok else 'BAD'}  {name:<20} expected -> {what}")
    bad = sum(not ok for _, ok, _ in results)
    print("SELF-TEST PASS" if not bad else f"SELF-TEST FAIL ({bad})")
    return bad


def intermediates(work: Path, root: Path):
    print(f"INTERMEDIATES: {work}  vs  {root}\n")
    n_fail = 0
    for wsub, rsub, names, first_col in INTERMEDIATE_MAP:
        wd, rd = work / wsub, root / rsub
        if not wd.exists():
            print(f"  (not built) {wsub}")
            continue
        for n in names or sorted(p.name for p in wd.glob("*.csv")):
            if not (rd / n).exists():
                print(f"  ----  {wsub}/{n}: no research copy to compare")
                continue
            if first_col is None:
                r = compare(wd / n, rd / n)
            else:
                a, b = (pd.read_csv(x, **READ_KW) for x in (wd / n, rd / n))
                keep = lambda d: [c for c in d.columns  # noqa: E731
                                  if not c.replace(".", "").isdigit()
                                  or float(c) >= first_col]
                r = compare_frames(a[keep(a)], b[keep(b)])
                r["note"] = (f"names + columns >= {first_col} only: "
                             + ("every value identical" if r["values"]
                                else r["note"] or r["where"]))
            ok = r["sha"] or r["values"]
            n_fail += not ok
            md = "" if r["max_abs"] is None else f"{r['max_abs']:.3e}"
            print(f"  {'PASS' if ok else 'FAIL'}  {wsub}/{n}  sha={'yes' if r['sha'] else 'no'}"
                  f" values={'yes' if r['values'] else 'no'} {md} {r['note'] or r['where']}")
    print(f"\n{n_fail} intermediate FAIL(s)")
    return n_fail


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--candidate-dir", type=Path)
    ap.add_argument("--reference-dir", type=Path, default=REFERENCE_DIR)
    ap.add_argument("--files", nargs="+", default=None)
    ap.add_argument("--report", type=Path, default=None, help="write the table as CSV")
    ap.add_argument("--allow-missing", action="store_true",
                    help="gate only the expected files that exist (a partial build); "
                         "list the rest as ABSENT")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--intermediates", type=Path, default=None,
                    help="a build_forcing.py --work-dir")
    ap.add_argument("--fairtofredi-root", type=Path, default=None)
    args = ap.parse_args()
    if args.self_test:
        sys.exit(1 if self_test() else 0)
    if args.intermediates:
        if not args.fairtofredi_root:
            sys.exit("ERROR: --intermediates needs --fairtofredi-root")
        sys.exit(1 if intermediates(args.intermediates,
                                    args.fairtofredi_root.expanduser()) else 0)
    if not args.candidate_dir:
        sys.exit("ERROR: --candidate-dir is required")
    sys.exit(1 if gate(args.candidate_dir, args.reference_dir, args.files,
                       args.report, args.allow_missing) else 0)


if __name__ == "__main__":
    main()
