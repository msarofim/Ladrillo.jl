#!/usr/bin/env python3
"""
check_targets_identity.py -- the gate: is a rebuild of the calibration targets
and drivers byte-identical to what Ladrillo ships (and, for the files it does
not ship, to the research repository's)?

    python tools/targets/check_targets_identity.py --work-dir /tmp/targets_build

Three groups, each compared by sha256:

  shipped        the eight data/ files (<work>/ladrillo/<path>) against
                 data/MANIFEST.csv's sha256, and their research-named twins
                 (<work>/<research path>) against MANIFEST's source_sha256.
                 Each reference data/ file is first checked against MANIFEST,
                 so a reference edited in place cannot pass.
  intermediate   four research files Ladrillo does not ship, against
                 tools/targets/reference.sha256.
  extraction     (only after build_targets.py --from-raw) the seven re-extracted
                 files (five shipped in data/sources/, two fetch-yourself)
                 against tools/targets/inputs.sha256.

There is NO tolerance. A file passes only on equal sha256, with ONE named
exception: the `commit` cell of gis_amp_shape_meta.csv records the commit of
the code that built it (d59ce12, research, in the shipped file; the Ladrillo
commit in a rebuild). For that file the gate requires every other byte to be
identical and reports the two commit labels (verdict PASS-LABEL). Nothing else
is exempt. A missing file FAILS.

On a failure it also reads both files with float_precision="round_trip" and
reports whether the values are equal (a formatting-only difference still
FAILS) or the largest |difference| and where.

--self-test mutation-tests this gate: a one-ulp change, a dropped row, a CRLF
line ending and a changed meta cell other than `commit` must FAIL; a changed
`commit` cell must PASS-LABEL; an untouched copy must PASS.
"""
from __future__ import annotations

import argparse
import csv
import io
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import target_inputs as TI  # noqa: E402
from build_targets import SHIPPED, INTERMEDIATES, LADRILLO_OUT, EXTRACTED_FILES  # noqa: E402

REPO = TI.REPO_ROOT
MANIFEST = REPO / "data" / "MANIFEST.csv"
REFERENCE = HERE / "reference.sha256"
LABEL_CELLS = {"constants/gis_amp_shape_meta.csv": "commit",
               "outputs/gis_amp_shape_meta.csv": "commit"}


def manifest():
    return {r["file"]: r for r in csv.DictReader(open(MANIFEST))}


def reference():
    out = {}
    for line in REFERENCE.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            sha, rel, *_ = line.split()
            out[rel] = sha
    return out


def value_diff(cand: Path, ref: Path) -> str:
    import numpy as np
    import pandas as pd
    try:
        a = pd.read_csv(cand, float_precision="round_trip")
        b = pd.read_csv(ref, float_precision="round_trip")
    except Exception as e:  # noqa: BLE001
        return f"unreadable as CSV ({type(e).__name__})"
    if list(a.columns) != list(b.columns):
        return f"columns differ: {list(a.columns)[:6]}... vs {list(b.columns)[:6]}..."
    if a.shape != b.shape:
        return f"shape {a.shape} vs {b.shape}"
    worst, where = 0.0, None
    for c in a.columns:
        x, y = a[c], b[c]
        if x.dtype.kind in "fi" and y.dtype.kind in "fi":
            d = (x.astype(float) - y.astype(float)).abs()
            both_nan = x.isna() & y.isna()
            d = d.where(~both_nan, 0.0)
            if d.isna().any():
                return f"NaN pattern differs in column {c}"
            if d.max() > worst:
                worst, where = float(d.max()), (c, int(d.idxmax()))
        elif not x.astype(str).equals(y.astype(str)):
            return f"text column {c} differs"
    if worst == 0.0:
        return "values equal (formatting-only difference)"
    return f"max |diff| {worst:.3e} in column {where[0]}, row {where[1]}"


def label_only_equal(cand: Path, ref: Path, col: str):
    """True + (cand label, ref label) if the two CSVs are byte-identical except the cells of column `col`."""
    ta, tb = cand.read_bytes(), ref.read_bytes()
    if b"\r" in ta or b"\r" in tb:
        return False, None
    ra = list(csv.reader(io.StringIO(ta.decode())))
    rb = list(csv.reader(io.StringIO(tb.decode())))
    if len(ra) != len(rb) or not ra or ra[0] != rb[0] or col not in ra[0]:
        return False, None
    k = ra[0].index(col)
    for x, y in zip(ra[1:], rb[1:]):
        if len(x) != len(y) or x[:k] + x[k + 1:] != y[:k] + y[k + 1:]:
            return False, None
    # re-serialise the candidate with the reference's label and require byte identity of the whole file
    swapped = ta.decode().splitlines(keepends=True)
    fixed = [swapped[0]]
    for line, x, y in zip(swapped[1:], ra[1:], rb[1:]):
        fixed.append(line.replace(x[k], y[k], 1) if x[k] != y[k] else line)
    if "".join(fixed).encode() != tb:
        return False, None
    return True, (ra[1][k], rb[1][k])


def compare(cand: Path, ref: Path, want_sha: str | None, label_col: str | None):
    """(verdict, detail). want_sha: the sha256 the reference must have (None: take the reference as is)."""
    if want_sha is not None:
        if not ref.exists():
            return "FAIL", f"reference {ref} missing"
        if TI.sha256(ref) != want_sha:
            return "FAIL", "REFERENCE does not match its recorded sha256 (edited in place?)"
    if not cand.exists():
        return "FAIL", "not built"
    got = TI.sha256(cand)
    ref_sha = want_sha if want_sha is not None else TI.sha256(ref)
    if got == ref_sha:
        return "PASS", "sha256 identical"
    if label_col and ref.exists():
        ok, labels = label_only_equal(cand, ref, label_col)
        if ok:
            return "PASS-LABEL", f"identical except the `{label_col}` label: {labels[0]} (rebuild) vs {labels[1]} (shipped)"
    detail = value_diff(cand, ref) if ref.exists() else "no reference file to diff"
    return "FAIL", detail


def gate(work: Path, research_root: Path | None, quiet=False):
    import json
    man, refs, pins = manifest(), reference(), TI.read_pins()
    prov = work / "BUILD_PROVENANCE.json"
    from_raw = bool(json.loads(prov.read_text()).get("from_raw")) if prov.exists() else False
    rows = []
    for dest, research_rel, _ in SHIPPED:
        m = man[dest]
        v, d = compare(work / LADRILLO_OUT / dest, REPO / "data" / dest, m["sha256"], LABEL_CELLS.get(dest))
        rows.append(("shipped", f"data/{dest}", v, d))
        # research-named twin vs MANIFEST source_sha256; diff against the research file when one is at hand
        rref = (research_root / research_rel) if research_root else (REPO / "data" / dest if m["transform"] == "copy"
                                                                       else Path("/nonexistent"))
        cand = work / research_rel
        if rref.exists() and TI.sha256(rref) == m["source_sha256"]:
            v, d = compare(cand, rref, m["source_sha256"], LABEL_CELLS.get(research_rel))
        else:
            got = TI.sha256(cand) if cand.exists() else None
            v, d = (("PASS", "sha256 identical") if got == m["source_sha256"] else
                    ("FAIL", "not built" if got is None else "sha256 differs (no research file at hand to diff)"))
            if v == "FAIL" and m["transform"] == "copy" and research_rel in LABEL_CELLS:
                v, d = compare(cand, REPO / "data" / dest, m["sha256"], LABEL_CELLS[research_rel])
        rows.append(("research", research_rel, v, d))
    for rel in INTERMEDIATES:
        cand = work / rel
        got = TI.sha256(cand) if cand.exists() else None
        if got == refs[rel]:
            v, d = "PASS", "sha256 identical"
        elif research_root and (research_root / rel).exists():
            v, d = compare(cand, research_root / rel, refs[rel], None)
        else:
            v, d = "FAIL", ("not built" if got is None else "sha256 differs from reference.sha256")
        rows.append(("intermediate", rel, v, d))
    if from_raw:
        for rel in EXTRACTED_FILES:
            cand = work / "sources" / rel
            if "sources/" + rel in pins:
                v, d = compare(cand, TI.SOURCES_DIR / rel, pins["sources/" + rel], None)
                rows.append(("extraction", f"data/sources/{rel}", v, d))
            else:
                got = TI.sha256(cand) if cand.exists() else None
                v, d = (("PASS", "sha256 identical") if got == pins["fetch/" + rel] else
                        ("FAIL", "not built" if got is None else "sha256 differs from inputs.sha256"))
                rows.append(("extraction", f"<inputs-dir>/{rel}", v, d))
    if not quiet:
        w = max(len(r[1]) for r in rows)
        for g, f, v, d in rows:
            print(f"  {g:12s} {f:{w}s}  {v:10s} {d}")
    n_fail = sum(r[2] == "FAIL" for r in rows)
    return rows, n_fail


def self_test():
    """Mutation-test the comparison itself on copies of shipped files."""
    tmp = Path(tempfile.mkdtemp(prefix="targets_gate_selftest_"))
    try:
        tgt = REPO / "data/observations/calibration_targets.csv"
        meta = REPO / "data/constants/gis_amp_shape_meta.csv"
        cases = []

        def case(name, ref, mutate, label_col, expect):
            c = tmp / f"{name}.csv"
            c.write_bytes(mutate(ref.read_bytes()))
            v, d = compare(c, ref, TI.sha256(ref), label_col)
            ok = v == expect
            cases.append(ok)
            print(f"  self-test {name:24s} expect {expect:10s} got {v:10s} {'ok' if ok else 'WRONG'}  ({d})")

        def one_ulp(b):
            import numpy as np
            lines = b.decode().split("\n")
            hdr = lines[0].split(",")
            j = hdr.index("dang")
            row = lines[1].split(",")
            x = float(row[j])
            row[j] = repr(float(np.nextafter(x, np.inf)))
            lines[1] = ",".join(row)
            return "\n".join(lines).encode()

        case("untouched", tgt, lambda b: b, None, "PASS")
        case("one_ulp_dang_1900", tgt, one_ulp, None, "FAIL")
        case("dropped_last_row", tgt, lambda b: b"\n".join(b.split(b"\n")[:-2]) + b"\n", None, "FAIL")
        case("crlf", tgt, lambda b: b.replace(b"\n", b"\r\n"), None, "FAIL")
        case("meta_commit_label", meta, lambda b: b.replace(b",d59ce12", b",abc1234"), "commit", "PASS-LABEL")
        case("meta_anchor_changed", meta, lambda b: b.replace(b"0.939826733054467", b"0.939826733054468"),
             "commit", "FAIL")
        case("meta_label_and_value", meta, lambda b: b.replace(b",d59ce12", b",abc1234").replace(
            b",40,south", b",39,south"), "commit", "FAIL")
        readme_ok = TI.readme_table_is_current()
        print(f"  self-test README sources table      {'matches target_inputs.py' if readme_ok else 'STALE'}")
        return all(cases) and readme_ok
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work-dir", help="a build_targets.py work directory")
    ap.add_argument("--research-root", help="SLR-RFF-BRICK checkout, to diff research-named files on failure")
    ap.add_argument("--report", help="write the table as CSV")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        ok = self_test()
        print("SELF-TEST " + ("PASS" if ok else "FAIL"))
        sys.exit(0 if ok else 1)
    if not a.work_dir:
        ap.error("--work-dir is required (or --self-test)")
    rows, n_fail = gate(Path(a.work_dir).resolve(), Path(a.research_root).resolve() if a.research_root else None)
    if a.report:
        with open(a.report, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["group", "file", "verdict", "detail"])
            w.writerows(rows)
    n = len(rows)
    print(f"TARGETS IDENTITY: {n - n_fail}/{n} pass" + (f", {n_fail} FAIL" if n_fail else "") +
          f" ({sum(r[2] == 'PASS-LABEL' for r in rows)} on the named provenance label)")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
