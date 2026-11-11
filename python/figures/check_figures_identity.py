#!/usr/bin/env python3
"""
check_figures_identity.py -- identity gate for the paper's Figs. 1-6: the outputs of python/figures/ against the
research renders that are in the paper.

  python python/figures/check_figures_identity.py <SLR-RFF-BRICK repo> [--pkgdir=outputs]

THE RULES
- PNG: byte-identical to the research file in <research>/figures/paper/. That holds only on the pinned runtime
  (requirements.txt: the matplotlib, NumPy and pandas versions the research renders were made with); elsewhere fonts and
  anti-aliasing move pixels. A PNG that is not byte-identical is reported with the number of differing pixels and their
  bounding box, so a rendering difference can be told from a data difference, and it FAILS either way.
- Caption sidecar (.caption.txt): identical text, except a difference declared in CAPTION_DIFF (a clause that names a
  research-only commit or file). Every declared difference must still occur, so a stale declaration fails too.
- CSV tables the figures read or write: same columns in the same order, same number of rows, and every cell EQUAL
  (numbers `==` as parsed with float_precision="round_trip", NaN == NaN; text exactly). The `provenance` column is
  excluded, as it names the producing script by construction.
Every difference is printed; the exit status is the number of failed checks.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GATE = "FIGURES-IDENTITY"
EXCLUDE_COLS = ["provenance"]

## (paper figure, package file under <pkgdir>/figures, research file under <research>/figures/paper)
FIGURES = [
    ("Fig. 1", "hindcast_components_L27.png", "hindcast_components_L27.png"),
    ("Fig. 2", "model_comparison_components_vv_L27_2100.png", "model_comparison_components_vv_L27_2100.png"),
    ("Fig. 2/3 sibling (2150)", "model_comparison_components_vv_L27_2150.png",
     "model_comparison_components_vv_L27_2150.png"),
    ("Fig. 3", "model_comparison_components_vv_L27_2300.png", "model_comparison_components_vv_L27_2300.png"),
    ("Fig. 4", "future_components_vv_L27_joint.png", "future_components_vv_L27_joint.png"),
    ("Fig. 5", "vv_gsic_ladrillo_L27_2300_joint.png", "vv_gsic_ladrillo_L27_2300_joint.png"),
    ("Fig. 6", "vv_responsiveness_L27.png", "vv_responsiveness_L27.png"),
]
## (what, package file under <pkgdir>/tables, research file under <research>/outputs)
TABLES = [
    ("Figs. 2, 3, 6 input: van Vuuren model comparison", "vv_model_comparison.csv", "vv_model_comparison_L27.csv"),
    ("Fig. 6 output: responsiveness", "vv_responsiveness.csv", "vv_responsiveness_L27.csv"),
]
## Declared caption differences: package file -> list of (research text, package text).
CAPTION_DIFF = {
    ## the research render git-logs its mean-GMST drivers; this repository has no research history, and its test suite
    ## checks the same files against data/MANIFEST.csv instead
    "vv_gsic_ladrillo_L27_2300_joint.png": [
        ("(driver commit 4854e59)",
         "(the forcing files are checked against data/MANIFEST.csv by the package's tests)")],
}


def check_png(pkg, ref):
    a, b = open(pkg, "rb").read(), open(ref, "rb").read()
    if a == b:
        return True, "byte-identical (%d bytes)" % len(a)
    try:
        from PIL import Image
        A = np.asarray(Image.open(pkg).convert("RGBA")).astype(int)
        B = np.asarray(Image.open(ref).convert("RGBA")).astype(int)
    except Exception as e:                                   # no Pillow: report bytes only
        return False, "bytes differ (%d vs %d); pixel comparison unavailable (%s)" % (len(a), len(b), e)
    if A.shape != B.shape:
        return False, "bytes differ; image size %s vs %s" % (A.shape[:2], B.shape[:2])
    d = np.any(A != B, axis=2)
    if not d.any():
        return False, "bytes differ but every pixel is equal (PNG metadata or compression)"
    ys, xs = np.nonzero(d)
    return False, ("%d of %d pixels differ, bounding box x %d-%d, y %d-%d"
                   % (d.sum(), d.size, xs.min(), xs.max(), ys.min(), ys.max()))


def check_caption(pkg_png, ref_png):
    pc, rc = os.path.splitext(pkg_png)[0] + ".caption.txt", os.path.splitext(ref_png)[0] + ".caption.txt"
    if not os.path.exists(pc) or not os.path.exists(rc):
        return False, "caption sidecar missing: %s" % (pc if not os.path.exists(pc) else rc)
    p, r = open(pc).read(), open(rc).read()
    for old, new in CAPTION_DIFF.get(os.path.basename(pkg_png), []):
        if old not in r or new not in p:
            return False, "declared caption difference no longer occurs: %r -> %r" % (old, new)
        r = r.replace(old, new)
    if p == r:
        n = len(CAPTION_DIFF.get(os.path.basename(pkg_png), []))
        return True, "identical" + (" after %d declared difference(s)" % n if n else "")
    return False, "caption differs:\n      package : %r\n      research: %r" % (p, r)


def check_csv(pkg, ref):
    P = pd.read_csv(pkg, float_precision="round_trip")
    R = pd.read_csv(ref, float_precision="round_trip")
    pc = [c for c in P.columns if c not in EXCLUDE_COLS]
    rc = [c for c in R.columns if c not in EXCLUDE_COLS]
    if pc != rc:
        return False, "columns differ:\n      package : %s\n      research: %s" % (pc, rc)
    if len(P) != len(R):
        return False, "rows differ: %d vs %d" % (len(P), len(R))
    bad = []
    for c in pc:
        a, b = P[c], R[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            ne = ~((a == b) | (a.isna() & b.isna()))
        else:
            ne = a.astype(str) != b.astype(str)
        if ne.any():
            bad.append("%s: %d cell(s), first at row %d (%r vs %r)"
                       % (c, int(ne.sum()), int(np.argmax(ne.values)), a[ne].iloc[0], b[ne].iloc[0]))
    if bad:
        return False, "cells differ:\n      " + "\n      ".join(bad)
    return True, "%d rows x %d columns equal" % (len(P), len(pc))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("research", help="the SLR-RFF-BRICK repository")
    ap.add_argument("--pkgdir", default=os.path.join(REPO, "outputs"))
    a = ap.parse_args()
    fails = 0
    print("[%s] package %s vs research %s" % (GATE, a.pkgdir, a.research))
    for what, pkg, ref in FIGURES:
        p = os.path.join(a.pkgdir, "figures", pkg)
        r = os.path.join(a.research, "figures", "paper", ref)
        for kind, fn in (("png", check_png), ("caption", check_caption)):
            if not os.path.exists(p) or not os.path.exists(r):
                ok, msg = False, "missing: %s" % (p if not os.path.exists(p) else r)
            else:
                ok, msg = fn(p, r)
            fails += not ok
            print("  %-4s %-24s %-8s %s" % ("PASS" if ok else "FAIL", what, kind, msg))
    for what, pkg, ref in TABLES:
        p, r = os.path.join(a.pkgdir, "tables", pkg), os.path.join(a.research, "outputs", ref)
        ok, msg = (check_csv(p, r) if os.path.exists(p) and os.path.exists(r)
                   else (False, "missing: %s" % (p if not os.path.exists(p) else r)))
        fails += not ok
        print("  %-4s %-24s %-8s %s" % ("PASS" if ok else "FAIL", what[:24], "csv", msg))
    print("[%s] %s: %d failed" % (GATE, "PASS" if fails == 0 else "FAIL", fails))
    sys.exit(fails)


if __name__ == "__main__":
    main()
