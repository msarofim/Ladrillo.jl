#!/usr/bin/env python3
"""Copy the run-control and extraction scripts behind data/comparison/ into comparison/, unchanged, and record each
file's source, source commit and sha256 in comparison/MANIFEST.csv (2026-11-22, for GMD's code and data policy:
"preprocessing, run control and postprocessing scripts covering every data processing action").

  python tools/copy_comparison_scripts.py --facts=<FACTS checkout> --research=<SLR-RFF-BRICK checkout>
  python tools/copy_comparison_scripts.py --check          # sha256 of every shipped file against the manifest

Files are read from git objects at the commit that produced the comparison data (`git show <commit>:<path>`), so a later
edit in either checkout cannot leak in. The one exception is the FACTS receipt arm's config, which was never committed
in the FACTS checkout; it is read from the working tree and its source commit is recorded as "untracked".
No RNG.
"""
import argparse, csv, hashlib, os, subprocess, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = os.path.join(REPO, "comparison")
MANIFEST = os.path.join(DEST, "MANIFEST.csv")

## FACTS: model code at tag v1.1.4 of https://github.com/radical-collaboration/facts; our additions are commits on top
## of it in a local checkout. d52b9ed4 (2026-10-08 18:21 -0400) is the state that wrote the FACTS file in
## data/comparison (SLR-RFF-BRICK e567ea7, 18:26 the same day).
FACTS_COMMIT = "d52b9ed4"
VV = ["vvVL", "vvLN", "vvL", "vvML", "vvM", "vvHL", "vvH"]
SSP = ["ssp126", "ssp245", "ssp585"]
FACTS_FILES = (["docker/Dockerfile", "build_shared_climate_nc.py", "build_shared_configs.py",
                "extract_facts_shared_components.py", "run_vv_cmip7_rerun_20261007.sh",
                "run_vv_zenodo_rerun_20261008.sh", "run_vv_facts_2300.sh", "run_ssp_facts_2300.sh"]
               + [f"experiments/global.shared.{s}.n200/config.yml" for s in VV + SSP]
               + [f"experiments/global.shared.{s}2300.n200/config.yml" for s in VV + SSP])
FACTS_UNTRACKED = ["experiments/global.shared.vvHL2300.n200.crate0/config.yml"]
## Research repository: the MAGICC-SLR extraction (commit that wrote the MAGICC files in data/comparison) and the
## receipt for the FACTS Greenland post-2100 numbers in Sect. 4 (46 and 112 cm at 2300 under High-to-Low).
RESEARCH_FILES = [("f7942da", "python/extract_magicc_components.py", "magicc/extract_magicc_components.py"),  # imported
                  ("f7942da", "python/extract_magicc_vv_components.py", "magicc/extract_magicc_vv_components.py"),
                  ("f7942da", "python/extract_magicc_vv_gmst.py", "magicc/extract_magicc_vv_gmst.py"),
                  ("e567ea7", "python/diag_facts_gis_extrap_receipt.py", "facts/receipt/diag_facts_gis_extrap_receipt.py"),
                  ("e567ea7", "outputs/diag_facts_gis_extrap_receipt_vvHL2300.csv",
                   "facts/receipt/diag_facts_gis_extrap_receipt_vvHL2300.csv")]
COLS = ["file", "source_repository", "source_commit", "source_path", "sha256"]


def sha(b):
    return hashlib.sha256(b).hexdigest()


def check():
    rows = list(csv.DictReader(open(MANIFEST)))
    bad = [r["file"] for r in rows if sha(open(os.path.join(DEST, r["file"]), "rb").read()) != r["sha256"]]
    listed = {r["file"] for r in rows}
    on_disk = {os.path.relpath(os.path.join(d, f), DEST) for d, _, fs in os.walk(DEST) for f in fs} - {"MANIFEST.csv",
                                                                                                   "README.md"}
    extra, missing = sorted(on_disk - listed), sorted(listed - on_disk)
    ok = not bad and not extra and not missing
    print(f"[COMPARISON-MANIFEST] {'PASS' if ok else 'FAIL'}: {len(rows)} files; sha256 mismatches {bad}; "
          f"unlisted {extra}; missing {missing}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--facts")
    ap.add_argument("--research")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        sys.exit(0 if check() else 1)
    rows = []

    def put(rel, data, repo, commit, src):
        out = os.path.join(DEST, rel)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        open(out, "wb").write(data)
        rows.append(dict(file=rel, source_repository=repo, source_commit=commit, source_path=src, sha256=sha(data)))

    show = lambda root, c, p: subprocess.run(["git", "-C", root, "show", f"{c}:{p}"], capture_output=True,
                                             check=True).stdout
    for p in FACTS_FILES:
        put("facts/" + p, show(a.facts, FACTS_COMMIT, p), "FACTS checkout (v1.1.4 + local commits)", FACTS_COMMIT, p)
    for p in FACTS_UNTRACKED:
        put("facts/" + p, open(os.path.join(a.facts, p), "rb").read(), "FACTS checkout (v1.1.4 + local commits)",
            "untracked", p)
    for c, p, rel in RESEARCH_FILES:
        put(rel, show(a.research, c, p), "SLR-RFF-BRICK", c, p)
    with open(MANIFEST, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["file"]))
    print(f"copied {len(rows)} files to comparison/")
    sys.exit(0 if check() else 1)


if __name__ == "__main__":
    main()
