"""Import Ladrillo v1.0 (L27) runtime data from the research repository, with provenance.

This is the ONE place the package's data/ directory was populated from SLR-RFF-BRICK. It is kept so that every file
in data/ can be traced to a source file and commit, and so the import can be re-run if the research repository
regenerates an input.

TRANSFORMS ARE LABEL-ONLY. The only edits are the glacier block renames
    SLOWP -> SLOWG   (slow-responding glacier regions: RGI 03, 09, 07, 06)
    FAST  -> FASTG   (the remaining 13 RGI regions)
applied to (a) column NAMES in a header line and (b) the `block` LABEL in glacier_block_constants.csv. Every numeric
field is copied byte-for-byte; data/MANIFEST.csv records the sha256 of source and destination and the transform, and
`verify()` re-checks that a renamed file differs from its source ONLY in the renamed tokens.

Usage:  python tools/import_l27_data.py --src ~/Documents/2026/CodeProjects/SLR-RFF-BRICK
"""
import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")

## Column-name renames. Anchored on a token boundary so FASTG cannot become FASTGG and words like FASTER, gis_fast,
## CUT_FASTDYN are untouched (the traps listed in the 2026-10-07 handoff, section 3F).
RENAME = [(re.compile(r"(?<![A-Za-z0-9])SLOWP(?![A-Za-z0-9])"), "SLOWG"),
          (re.compile(r"(?<![A-Za-z0-9])FAST(?![A-Za-z0-9])"), "FASTG")]
## The same tokens as a column-name SUFFIX (gic_a_FAST, gsic_slowp is a Mimi variable, not data).
RENAME_SUFFIX = [(re.compile(r"_SLOWP\b"), "_SLOWG"), (re.compile(r"_FAST\b"), "_FASTG")]

FORCING_KEYS = ["ssp126", "ssp245", "ssp585"] + ["vv" + m for m in ("VL", "LN", "L", "ML", "M", "HL", "H")]

## (destination under data/, source under the research repo, transform)
FILES = [
    ("posterior/ladrillo_L27_posterior.csv", "data/MimiBRICK/parameters_subsample_brick_mengel_L27.csv", "header"),
    ("drivers/glacier_block_temperature.csv", "data/observations/t_glac_blocks.csv", "header"),
    ("drivers/greenland_zone_temperature.csv", "data/observations/t_gis_zones.csv", "copy"),
    ("constants/glacier_block_constants.csv", "outputs/extc_block_constants.csv", "block_label"),
    ("constants/gis_amp_shape.csv", "outputs/gis_amp_shape.csv", "copy"),
    ("constants/gis_amp_shape_meta.csv", "outputs/gis_amp_shape_meta.csv", "copy"),
    ("constants/central_row.csv", "outputs/recalib_central_row.csv", "copy"),
    ("constants/paleo_fastdyn_draws.csv", "outputs/paleo_fastdyn_draws.csv", "copy"),
    ("constants/paleo_dais_marginals.csv", "outputs/paleo_dais_marginals.csv", "copy"),
    ("observations/calibration_targets.csv", "outputs/recalib_targets_ext.csv", "copy"),
    ("forcing/fair_mean_gmst_ssp245harm.csv", "data/observations/fair_mean_gmst_ssp245harm.csv", "copy"),
    ("forcing/fair_mean_ohc_ssp245harm.csv", "data/observations/fair_mean_ohc_ssp245harm.csv", "copy"),
    ("forcing/fair_cube_gmst_ssp585_spliced.csv", "data/observations/fair_cube_gmst_ssp585_spliced.csv", "copy"),
    ## BRICK 2.0 comparison arm: the post-PR#93 MimiBRICK v2.0.0 posterior subsample (10,000 draws; linear-units
    ## antarctic_precip0, see brick20.jl). Redistribution to be confirmed with its author before the repo is public.
    ("brick20/parameters_subsample_brick.csv", "data/MimiBRICK/parameters_subsample_brick.csv", "copy"),
] + [("forcing/fair_%s_%s_%s.csv" % (kind, var, k) if kind == "mean" else
      "forcing/fair_%s_%s_%s_raw.csv" % (kind, var, k),
      "data/observations/fair_%s_%s_%s.csv" % (kind, var, k) if kind == "mean" else
      "data/observations/fair_%s_%s_%s_raw.csv" % (kind, var, k), "copy")
     for k in FORCING_KEYS for kind in ("mean", "cube") for var in ("gmst", "ohc")]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rename_tokens(s):
    for rx, new in RENAME_SUFFIX + RENAME:
        s = rx.sub(new, s)
    return s


def transform(src, dst, how):
    with open(src, "rb") as f:
        raw = f.read()
    if how == "copy":
        out = raw
    elif how == "header":
        nl = raw.index(b"\n")
        out = rename_tokens(raw[:nl].decode("utf-8")).encode("utf-8") + raw[nl:]
    elif how == "block_label":
        lines = raw.split(b"\n")
        for i in range(1, len(lines)):
            if lines[i]:
                head, sep, rest = lines[i].partition(b",")
                lines[i] = rename_tokens(head.decode("utf-8")).encode("utf-8") + sep + rest
        out = b"\n".join(lines)
    else:
        raise SystemExit("unknown transform %s" % how)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "wb") as f:
        f.write(out)


def verify(src, dst, how):
    """A renamed file must differ from its source ONLY in the renamed tokens: undoing the renames on the destination
    must reproduce the source byte-for-byte."""
    a, b = open(src, "rb").read(), open(dst, "rb").read()
    if how == "copy":
        return a == b
    back = b.decode("utf-8").replace("SLOWG", "SLOWP")
    back = re.sub(r"(?<![A-Za-z0-9])FASTG(?![A-Za-z0-9])", "FAST", back)
    return back.encode("utf-8") == a and b"FASTGG" not in b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    a = ap.parse_args()
    src_root = os.path.expanduser(a.src)
    commit = subprocess.run(["git", "-C", src_root, "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    rows, bad = [], 0
    for dst_rel, src_rel, how in FILES:
        s, d = os.path.join(src_root, src_rel), os.path.join(DATA, dst_rel)
        if not os.path.exists(s):
            sys.exit("missing source %s" % s)
        transform(s, d, how)
        ok = verify(s, d, how)
        bad += not ok
        rows.append(dict(file=dst_rel, source=src_rel, source_commit=commit, transform=how,
                         source_sha256=sha(s), sha256=sha(d), label_only_check="PASS" if ok else "FAIL"))
        print("%-55s <- %-60s %-11s %s" % (dst_rel, src_rel, how, "PASS" if ok else "FAIL"))
    with open(os.path.join(DATA, "MANIFEST.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    if bad:
        sys.exit("%d file(s) failed the label-only check" % bad)
    print("imported %d files from SLR-RFF-BRICK @ %s; data/MANIFEST.csv written" % (len(rows), commit))


if __name__ == "__main__":
    main()
