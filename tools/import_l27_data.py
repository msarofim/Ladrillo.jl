"""Import Ladrillo's runtime data (calibration L27) from the research repository, with provenance.

This is the ONE place the package's data/ directory was populated from SLR-RFF-BRICK. It is kept so that every file
in data/ can be traced to a source file and commit, and so the import can be re-run if the research repository
regenerates an input.

TRANSFORMS ARE LABEL-ONLY. The only edits are the glacier block renames
    SLOWP -> SLOWG   (slow-responding glacier regions: RGI 03, 09, 07, 06)
    FAST  -> FASTG   (the remaining 13 RGI regions)
applied to (a) column NAMES in a header line, (b) the `block` LABEL in glacier_block_constants.csv and (c) for
transform "labels", every token in the file (ladrillo_priors_L27.csv carries the names in its `name` column and inside
its `prior_form` text). Every numeric field is copied byte-for-byte; data/MANIFEST.csv records the sha256 of source and
destination and the transform, and `verify()` re-checks that a renamed file differs from its source ONLY in the renamed
tokens.

Usage:  python tools/import_l27_data.py --src <SLR-RFF-BRICK checkout>
        python tools/import_l27_data.py --src ... --only constants/ladrillo_priors_L27.csv constants/paleo_geo_prior_ton.csv
--only imports just the named destinations: it refuses a source with uncommitted changes, rewrites (or appends) only
those files' MANIFEST rows, and keeps every other MANIFEST row byte-for-byte (its file was imported at its own commit).
"""
import argparse
import csv
import hashlib
import io
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")

## Column-name renames. Anchored on a token boundary so FASTG cannot become FASTGG and words like FASTER, gis_fast,
## CUT_FASTDYN are untouched.
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
    ## the glacier target with the RGI-region-19 seam removed (what the hindcast glacier scope is compared against)
    ("observations/calibration_targets_gsic_adj.csv", "outputs/recalib_targets_ext_gsicadj.csv", "copy"),
    ("forcing/fair_mean_gmst_ssp245harm.csv", "data/observations/fair_mean_gmst_ssp245harm.csv", "copy"),
    ("forcing/fair_mean_ohc_ssp245harm.csv", "data/observations/fair_mean_ohc_ssp245harm.csv", "copy"),
    ("forcing/fair_cube_gmst_ssp585_spliced.csv", "data/observations/fair_cube_gmst_ssp585_spliced.csv", "copy"),
    ## BRICK 2.0 comparison arm: the MimiBRICK v2.0.0 joint SNEASY-BRICK posterior subsample (Zenodo 20592337,
    ## 10,000 draws; log-space antarctic_precip0, see brick20.jl). Not shipped: tools/fetch_brick20_posterior.sh
    ## fetches it from Zenodo.
    ("brick20/parameters_subsample_sneasybrick.csv", "data/MimiBRICK/parameters_subsample_sneasybrick_zenodo20592337.csv",
     "copy"),
    ## Tables A1/A2 (python/tables/): the calibrator's own dump of the L27 priors (calibrate_mcmc_ext.jl --dump-priors,
    ## the production flags) and the DAIS paleo-ensemble joint prior of the seven geometry parameters
    ("constants/ladrillo_priors_L27.csv", "outputs/ladrillo_priors_L27.csv", "labels"),
    ("constants/paleo_geo_prior_ton.csv", "outputs/paleo_geo_prior_ton.csv", "copy"),
    ## The calibration (src/calibration.jl, scripts/calibrate.jl): the four over-dispersed start rows (one per chain,
    ## seeds 2026-2029), the proposal-covariance seed (adapted by the L26 run; it cannot be rebuilt without re-running
    ## L26) and the Greenland amplification prior table (calibration L27 reads the south/full row).
    ("calibration/overdispersed_starts_L27.csv", "outputs/mcmc/overdispersed_starts_L27.csv", "header"),
    ("calibration/adapted_cov_L26_named.csv", "outputs/mcmc/adapted_cov_L26_named.csv", "header"),
    ("calibration/gis_amp_prior.csv", "outputs/gis_amp_prior.csv", "copy"),
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
    elif how == "labels":
        out = rename_tokens(raw.decode("utf-8")).encode("utf-8")
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


## every token that carries a renamed block name, for checking the renames by eye
TOKEN = re.compile(r"[A-Za-z0-9_]*(?:SLOWP|SLOWG|FAST|FASTG)[A-Za-z0-9_]*")


def renamed_tokens(src, dst):
    """(source token, destination token, count) for every token the transform changed, in order of appearance."""
    a = TOKEN.findall(open(src, encoding="utf-8").read())
    b = TOKEN.findall(open(dst, encoding="utf-8").read())
    if len(a) != len(b):
        return [("<token count %d>" % len(a), "<token count %d>" % len(b), 1)]
    seen = {}
    for x, y in zip(a, b):
        if x != y:
            seen[(x, y)] = seen.get((x, y), 0) + 1
    return [(x, y, n) for (x, y), n in seen.items()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--only", nargs="+", metavar="DEST",
                    help="import only these destinations (paths under data/); keep every other MANIFEST row as is")
    a = ap.parse_args()
    src_root = os.path.expanduser(a.src)
    commit = subprocess.run(["git", "-C", src_root, "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    todo = FILES
    if a.only:
        unknown = sorted(set(a.only) - {d for d, _, _ in FILES})
        if unknown:
            sys.exit("--only: not in FILES: %s" % unknown)
        todo = [f for f in FILES if f[0] in a.only]
        ## the row records the source commit, so the source must BE that commit's file
        dirty = subprocess.run(["git", "-C", src_root, "status", "--porcelain", "--"] + [s for _, s, _ in todo],
                               capture_output=True, text=True).stdout.strip()
        if dirty:
            sys.exit("sources have uncommitted changes, refusing:\n%s" % dirty)
    rows, bad = [], 0
    for dst_rel, src_rel, how in todo:
        s, d = os.path.join(src_root, src_rel), os.path.join(DATA, dst_rel)
        if not os.path.exists(s):
            sys.exit("missing source %s" % s)
        transform(s, d, how)
        ok = verify(s, d, how)
        bad += not ok
        rows.append(dict(file=dst_rel, source=src_rel, source_commit=commit, transform=how,
                         source_sha256=sha(s), sha256=sha(d), label_only_check="PASS" if ok else "FAIL"))
        print("%-55s <- %-60s %-11s %s" % (dst_rel, src_rel, how, "PASS" if ok else "FAIL"))
        if a.only and how != "copy":
            for x, y, n in renamed_tokens(s, d):
                print("      renamed %-28s -> %-28s x%d" % (x, y, n))
    man = os.path.join(DATA, "MANIFEST.csv")
    if not a.only:
        with open(man, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    else:
        ## keep every other row's bytes; replace an imported file's row in place, append a new one at the end
        raw = open(man, "rb").read()
        eol = b"\r\n" if b"\r\n" in raw else b"\n"
        lines = raw.split(eol)
        assert lines[-1] == b"", "MANIFEST.csv does not end in a line break"
        lines = lines[:-1]
        fields = next(csv.reader([lines[0].decode("utf-8")]))
        assert fields == list(rows[0].keys()), "MANIFEST.csv columns changed: %s" % fields

        def enc(r):
            buf = io.StringIO()
            csv.DictWriter(buf, fieldnames=fields, lineterminator="").writerow(r)
            return buf.getvalue().encode("utf-8")
        where = {next(csv.reader([ln.decode("utf-8")]))[0]: i for i, ln in enumerate(lines) if i > 0}
        for r in rows:
            if r["file"] in where:
                lines[where[r["file"]]] = enc(r)
            else:
                lines.append(enc(r))
        with open(man, "wb") as f:
            f.write(eol.join(lines) + eol)
    if bad:
        sys.exit("%d file(s) failed the label-only check" % bad)
    print("imported %d files from SLR-RFF-BRICK @ %s; data/MANIFEST.csv %s" %
          (len(rows), commit, "updated (those rows only)" if a.only else "written"))


if __name__ == "__main__":
    main()
