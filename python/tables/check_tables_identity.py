#!/usr/bin/env python3
"""
check_tables_identity.py -- identity gate for Tables 4, A1 and A2: the outputs of python/tables/ against the research
outputs that produced the paper.

  python python/tables/check_tables_identity.py <SLR-RFF-BRICK repo> [--pkgdir=outputs/tables]
  python python/tables/check_tables_identity.py --refdir=test/reference [--pkgdir=DIR] [--rtol=0]

The first form compares all five outputs with the research repository's files. The second compares the same five
against the copies tools/build_test_reference.py makes in test/reference/; test/runtests.jl uses it. (The PCA CSV is
the one that sees a 1-ulp input change: Table A2's markdown prints 3 significant figures.)

THE RULES
- Markdown (.md): the text must be identical line by line, except
    (a) the trailing provenance comment (`<!-- ... -->`, the last non-empty line), which names the producing script and
        its inputs and so differs by construction; both files must have one, and both are printed;
    (b) a line that becomes identical once the research line's glacier-block labels are renamed (RENAME below).
- CSV: same columns in the same order and the same number of rows; then every cell must be EQUAL to the research
  cell after the research cell's glacier-block labels are renamed: numbers `==` as parsed by pandas with
  float_precision="round_trip" (NaN == NaN), text exactly. The `provenance` column is excluded (it differs by
  construction). A research cell holding a renamed label must appear RENAMED in the package (a package that kept the
  old name fails).
- Every difference is printed. There is no tolerance unless --rtol is given; test/runtests.jl passes one only off the
  pinned Python runtime, and the number of values that are not bit-identical is printed either way.
- No package output may contain FASTGG (a rename applied twice).

RENAME is the package's glacier-block rename, SLOWP -> SLOWG and FAST -> FASTG as whole tokens (an underscore counts
as a separator, so gic_a_FAST -> gic_a_FASTG; CUT_FASTDYN, FASTER and the Greenland gis_fast* / gis_alpha_f names are
untouched). It is written out here, not imported from tools/import_l27_data.py, so that a wrong rename in the import
cannot also be the gate's mapping. Every renamed research label is listed.
"""
import argparse
import math
import os
import re
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_PKGDIR = os.path.join(REPO, "outputs", "tables")
GATE = "TABLES-IDENTITY"
EXCLUDE_COLS = ["provenance"]
DOUBLE_RENAME = "FASTGG"
RENAME = {"SLOWP": "SLOWG", "FAST": "FASTG"}
_RENAME_RX = re.compile(r"(?<![A-Za-z0-9])(" + "|".join(RENAME) + r")(?![A-Za-z0-9])")
PROVENANCE_RX = re.compile(r"^<!--.*-->$")

## (output, kind, research file under the research repo, reference file under test/reference or None, package file)
PAIRS = [
    ("Table A2 basis: Antarctic-block PCA", "csv", "outputs/diag_ais_block_pca_L27_sub10k.csv",
     "research_diag_ais_block_pca_L27_sub10k.csv",
     "diag_ais_block_pca.csv"),
    ("Table A2", "md", "outputs/ladrillo_table_a2_L27_sub10k.md", "research_ladrillo_table_a2_L27_sub10k.md",
     "ladrillo_table_a2.md"),
    ("Table A1 (csv)", "csv", "outputs/ladrillo_prior_posterior_L27.csv", "research_ladrillo_prior_posterior_L27.csv",
     "ladrillo_prior_posterior.csv"),
    ("Table A1 (md)", "md", "outputs/ladrillo_prior_posterior_L27.md", "research_ladrillo_prior_posterior_L27.md",
     "ladrillo_prior_posterior.md"),
    ("Table 4: hindcast scorecard", "csv", "outputs/scope_ladrillo_vs_brick20_scorecard_L27.csv",
     "research_scope_ladrillo_vs_brick20_scorecard_L27.csv", "ladrillo_vs_brick20_scorecard.csv"),
]


def rename(s):
    return _RENAME_RX.sub(lambda m: RENAME[m.group(1)], s)


def isnan(x):
    return isinstance(x, float) and math.isnan(x)


def compare_md(res, pkg, renamed):
    """Return the list of differences (strings)."""
    a = open(res, encoding="utf-8").read().split("\n")
    b = open(pkg, encoding="utf-8").read().split("\n")
    diffs = []
    prov = []
    for name, lines in (("research", a), ("package", b)):
        k = max(i for i, ln in enumerate(lines) if ln.strip()) if any(ln.strip() for ln in lines) else -1
        if k < 0 or not PROVENANCE_RX.match(lines[k]):
            diffs.append(f"{name} file has no trailing provenance comment")
            prov.append(None)
        else:
            prov.append(lines.pop(k))
    print(f"      provenance (excluded)  research: {prov[0]}")
    print(f"      provenance (excluded)  package : {prov[1]}")
    if len(a) != len(b):
        diffs.append(f"line count {len(a)} (research) vs {len(b)} (package)")
    for i, (x, y) in enumerate(zip(a, b), start=1):
        if x == y:
            continue
        if rename(x) == y:
            renamed.append(f"line {i}")
            continue
        j = next((k for k, (p, q) in enumerate(zip(x, y)) if p != q), min(len(x), len(y)))
        diffs.append(f"line {i}, first difference at character {j + 1}: research ...{x[max(0, j - 30):j + 30]!r}... "
                     f"vs package ...{y[max(0, j - 30):j + 30]!r}...\n          research: {x}\n          package : {y}")
    return diffs


def compare_csv(res, pkg, rtol, renamed, stats):
    a = pd.read_csv(res, float_precision="round_trip")
    b = pd.read_csv(pkg, float_precision="round_trip")
    a = a.drop(columns=[c for c in EXCLUDE_COLS if c in a.columns])
    b = b.drop(columns=[c for c in EXCLUDE_COLS if c in b.columns])
    acols = [rename(c) for c in a.columns]
    if acols != list(b.columns):
        return [f"columns differ:\n          research (renamed): {acols}\n          package           : {list(b.columns)}"]
    if len(a) != len(b):
        return [f"row count {len(a)} (research) vs {len(b)} (package)"]
    key = a.columns[0]
    diffs = []
    for ca, cb in zip(a.columns, b.columns):
        for i, (x, y) in enumerate(zip(a[ca].tolist(), b[cb].tolist())):
            where = f"row {i + 1} ({key}={a[key].iloc[i]}), column {cb}"
            if isinstance(x, str) or isinstance(y, str):
                xm = rename(x) if isinstance(x, str) else x
                if isinstance(x, str) and xm != x:
                    renamed.append(f"{x} -> {xm}")
                if xm == y and type(xm) is type(y):
                    continue
                diffs.append(f"{where}: research {x!r} (renamed {xm!r}) vs package {y!r}")
                continue
            stats["n"] += 1
            if (isnan(x) and isnan(y)) or (x == y and not isnan(x)):
                stats["exact"] += 1
                continue
            num = isinstance(x, (int, float, np.integer, np.floating)) and isinstance(y, (int, float, np.integer, np.floating))
            rel = abs(x - y) / max(abs(x), abs(y)) if num and not (isnan(x) or isnan(y)) else math.inf
            stats["maxrel"] = max(stats["maxrel"], rel)
            if rtol > 0 and rel <= rtol:
                continue                                           # counted as not bit-identical, printed below
            diffs.append(f"{where}: research {x!r} vs package {y!r} (rel diff {rel:.3e})")
    return diffs


def main():
    ap = argparse.ArgumentParser(description="identity gate for Tables 4, A1 and A2 (python/tables/)")
    ap.add_argument("research", nargs="?", help="the SLR-RFF-BRICK repository")
    ap.add_argument("--refdir", help="compare against test/reference/ instead of the research repository")
    ap.add_argument("--pkgdir", default=DEFAULT_PKGDIR, help="the package outputs (default outputs/tables)")
    ap.add_argument("--rtol", type=float, default=0.0, help="relative tolerance on numbers (default 0: exact)")
    a = ap.parse_args()
    if bool(a.research) == bool(a.refdir):
        ap.error("give exactly one of <research repo> or --refdir")
    print(f"[{GATE}] python {sys.version.split()[0]} numpy {np.__version__} pandas {pd.__version__} | package outputs "
          f"{a.pkgdir} | against {a.research or a.refdir} | numbers {'EXACT' if a.rtol == 0 else f'to rtol {a.rtol:g}'}")
    bad = 0
    for label, kind, res_rel, ref_name, pkg_name in PAIRS:
        if a.refdir:
            if ref_name is None:
                print(f"  {label:40s} not in the reference set (compared only against the research repository)")
                continue
            res = os.path.join(a.refdir, ref_name)
        else:
            res = os.path.join(os.path.expanduser(a.research), res_rel)
        pkg = os.path.join(a.pkgdir, pkg_name)
        print(f"  {label:40s} {os.path.basename(res)}  vs  {pkg_name}")
        missing = [p for p in (res, pkg) if not os.path.isfile(p)]
        if missing:
            print(f"      MISSING: {missing}  FAIL")
            bad += 1
            continue
        renamed, stats = [], dict(n=0, exact=0, maxrel=0.0)
        diffs = compare_md(res, pkg, renamed) if kind == "md" else compare_csv(res, pkg, a.rtol, renamed, stats)
        if DOUBLE_RENAME in open(pkg, encoding="utf-8").read():
            diffs.append(f"package file contains {DOUBLE_RENAME}")
        for d in diffs:
            print(f"      DIFF {d}")
        if renamed:
            uniq = sorted(set(renamed)) if kind == "csv" else renamed
            print(f"      {len(renamed)} research {'cells' if kind == 'csv' else 'lines'} matched only after the rename: "
                  + ", ".join(uniq))
        if kind == "csv":
            print(f"      {stats['n']} numbers, {stats['exact']} bit-identical, {stats['n'] - stats['exact']} not, "
                  f"max rel diff {stats['maxrel']:.3e}")
        print(f"      {'IDENTICAL' if not diffs else f'FAIL ({len(diffs)} differences)'}")
        bad += len(diffs)
    print(f"[{GATE}] {'PASS' if bad == 0 else 'FAIL'}")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
