#!/usr/bin/env python3
"""
import_sources.py -- copy the shipped source files from the research repository
into data/sources/ (maintainer tool; the one place data/sources/ was populated).

Every file is copied byte-for-byte, with two exceptions, each a named transform:

  zip-member   glambie/calendar_years/*.csv are members of the research repo's
               glambie_data.zip (glambie_results_20240716/calendar_years/), read
               out of the zip unchanged.
  row-subset   gmip3/gmip3_block_ladder_L27.csv keeps the header and the rows of
               the research ladder cache whose `key` is one of the three region
               sets Ladrillo's glacier reservoirs use (R19, SLOWG, FASTG). The
               cache also holds rows for region sets tried and dropped during the
               2026-08 development; ladrillo_data.py reads only these three keys
               and fails loudly on a missing one. Every kept line is unchanged.

Only `sources/` (shipped) files are copied. Files with no research path (the
extractions grace_lws_domain_monthly.csv and greenland_zone_series_HadCRUT5.csv)
are written by tools/targets/extract/ instead.

Usage: python import_sources.py --research <SLR-RFF-BRICK checkout>
       python import_sources.py --write-pins --inputs-dir DIR
  --write-pins rewrites inputs.sha256: `sources/` lines from data/sources/,
  `fetch/` lines from DIR (laid out by fetch_inputs.py), `raw/` lines from DIR/raw/.
"""
import argparse
import os
import subprocess
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import target_inputs as TI  # noqa: E402

GLAMBIE_ZIP = "data/observations/raw/glambie_data.zip"
LADDER_KEYS = {"19", "03-06-07-09", "01-02-04-08-10-11-12-13-14-15-16-17-18"}   # R19, SLOWG, FASTG


def copy_one(research_root, rel, research):
    dst = TI.SOURCES_DIR / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if research.startswith("glambie_data.zip:"):
        member = research.split(":", 1)[1]
        with zipfile.ZipFile(os.path.join(research_root, GLAMBIE_ZIP)) as z:
            raw = z.read(member)
        how = "zip-member"
    elif rel == "gmip3/gmip3_block_ladder_L27.csv":
        lines = open(os.path.join(research_root, research), "rb").read().split(b"\n")
        assert lines[-1] == b"", "ladder cache does not end in a newline"
        keep = [lines[0]] + [ln for ln in lines[1:-1] if ln.split(b",", 1)[0].decode() in LADDER_KEYS]
        assert len(keep) == 1 + 4 * len(LADDER_KEYS), f"expected 4 rungs x 3 keys, got {len(keep) - 1} rows"
        raw = b"\n".join(keep) + b"\n"
        how = "row-subset"
    else:
        raw = open(os.path.join(research_root, research), "rb").read()
        how = "copy"
    dst.write_bytes(raw)
    return how


def write_pins(inputs_dir):
    """Rewrite inputs.sha256 from data/sources/ (sources/), DIR (fetch/) and DIR/raw (raw/)."""
    rows = []
    shipped = {rel for _, rel, _, _ in TI.files("sources")}
    on_disk = {p.relative_to(TI.SOURCES_DIR).as_posix() for p in TI.SOURCES_DIR.rglob("*")
               if p.is_file() and p.name != "README.md"}
    assert on_disk == shipped, f"data/sources/ differs from the table: extra {on_disk - shipped}, missing {shipped - on_disk}"
    for rel in sorted(shipped):
        rows.append(f"{TI.sha256(TI.SOURCES_DIR / rel)}  sources/{rel}")
    for _, rel, _, _ in sorted(TI.files("fetch"), key=lambda t: t[1]):
        rows.append(f"{TI.sha256(os.path.join(inputs_dir, rel))}  fetch/{rel}")
    for r in TI.RAW:
        rows.append(f"{TI.sha256(os.path.join(inputs_dir, 'raw', r['rel']))}  raw/{r['rel']}")
    head = ["# sha256 of every input of tools/targets/ (target_inputs.py says what each is).",
            "#   sources/<rel>  data/sources/<rel>             shipped",
            "#   fetch/<rel>    <--inputs-dir>/<rel>           fetch yourself (fetch_inputs.py); not shipped",
            "#   raw/<rel>      <--inputs-dir>/raw/<rel>       raw products, read only by the maintainer extractions"]
    TI.SHA_FILE.write_text("\n".join(head + rows) + "\n")
    print(f"wrote {TI.SHA_FILE} ({len(rows)} lines)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--research")
    ap.add_argument("--write-pins", action="store_true")
    ap.add_argument("--inputs-dir", help="with --write-pins: the fetch_inputs.py directory")
    a = ap.parse_args()
    if a.write_pins:
        if not a.inputs_dir:
            ap.error("--write-pins needs --inputs-dir")
        write_pins(a.inputs_dir)
        return
    if not a.research:
        ap.error("--research is required")
    root = os.path.expanduser(a.research)
    commit = subprocess.run(["git", "-C", root, "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    for _, rel, research, key in TI.files("sources"):
        if research is None:
            continue
        how = copy_one(root, rel, research)
        last = subprocess.run(["git", "-C", root, "log", "-1", "--format=%h", "--",
                               research.split(":")[0] if ":" not in research else GLAMBIE_ZIP],
                              capture_output=True, text=True).stdout.strip()
        print(f"{rel:68s} <- {research[:70]:70s} {how:10s} last commit {last}")
    print(f"research HEAD {commit}")


if __name__ == "__main__":
    main()
