"""
fetch_inputs.py -- download the PUBLIC inputs of the forcing build into an
--inputs-dir laid out as forcing_inputs.INPUTS says, and verify every file
against inputs.sha256.

    python fetch_inputs.py --inputs-dir INPUTS              # PLAN ONLY (default)
    python fetch_inputs.py --inputs-dir INPUTS --download   # actually fetch
    python fetch_inputs.py --inputs-dir INPUTS --verify     # hash what is there

Without --download nothing touches the network: the script prints what it
would fetch, from where, and to which path.

WHAT IT CANNOT FETCH. Two inputs are not public (forcing_inputs.INPUTS
`public=False`):
  * the IIASA scenariomip-2026-prerelease pull (gated database), which feeds the
    2024-2100 block of every van Vuuren marker; and
  * EDF / Lisa Rennels' AR6 SSP2-4.5 1750-2300 file, whose post-2100 RATIOS
    shape the ssp245harm calibration forcing's 2101-2300 tail.
The script names them and the path each must be placed at; it does not invent a
source. Until they are placed, build_forcing.py refuses to run (its input check
fails), rather than building something that is not the shipped forcing.

RCMIP v5.1.0 is fetched with fair's own pooch.retrieve(url, known_hash) call, so
it lands in the cache fair.fill_from_rcmip() reads; a copy is kept under
<inputs-dir>/rcmip/ so `shasum -a 256 -c inputs.sha256` works from INPUTS.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import forcing_inputs as fi  # noqa: E402

USER_AGENT = "Ladrillo-forcing-fetch/1.0"


def _download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh)
    tmp.replace(dest)


def fetch_one(key, inputs_dir):
    spec = fi.INPUTS[key]
    dest = Path(inputs_dir) / spec["rel"]
    if key in fi.RCMIP_FILES:
        import pooch
        fn, known = fi.RCMIP_FILES[key]
        cached = Path(pooch.retrieve(url=fi.RCMIP_BASE + fn, known_hash=known))
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cached, dest)
        return dest
    if "zip_member_suffix" in spec:
        with tempfile.TemporaryDirectory() as td:
            z = Path(td) / "record.zip"
            _download(spec["url"], z)
            with zipfile.ZipFile(z) as zf:
                hits = [n for n in zf.namelist()
                        if n.endswith(spec["zip_member_suffix"])]
                if len(hits) != 1:
                    sys.exit(f"ERROR: {key}: expected one member ending "
                             f"{spec['zip_member_suffix']!r}, found {hits}")
                dest.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(hits[0]) as src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out)
        return dest
    _download(spec["url"], dest)
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs-dir", type=Path, required=True)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--download", action="store_true",
                   help="actually fetch the public inputs (network)")
    g.add_argument("--verify", action="store_true",
                   help="only hash what is already in --inputs-dir")
    args = ap.parse_args()
    want = fi.expected_sha256()

    if args.verify:
        fails = fi.verify(inputs_dir=args.inputs_dir)
        print(f"\n{len(fi.INPUTS) - len(fails)} of {len(fi.INPUTS)} inputs verified")
        sys.exit(1 if fails else 0)

    print(f"{'DOWNLOADING' if args.download else 'PLAN (no network; pass --download)'}"
          f" into {args.inputs_dir}\n")
    private = []
    for key, spec in fi.INPUTS.items():
        dest = Path(args.inputs_dir) / spec["rel"]
        if not spec["public"]:
            private.append((key, spec, dest))
            continue
        print(f"  {key:<22} {spec['url']}\n  {'':<22} -> {dest}")
        if args.download:
            got = fetch_one(key, args.inputs_dir)
            h = fi.sha256(got)
            ok = h == want[spec["rel"]]
            print(f"  {'':<22} sha256 {'PASS' if ok else 'FAIL ' + h}")
            if not ok:
                print(f"  {'':<22} ⚠ the fetched file is NOT the one the shipped "
                      f"forcing was built from; build_forcing.py will refuse it.")
    print("\nNOT PUBLIC -- obtain these yourself and place them at the path shown:")
    for key, spec, dest in private:
        print(f"  {key:<22} {spec['source']}\n  {'':<22} -> {dest}\n"
              f"  {'':<22} sha256 {want[spec['rel']]}")
    if args.download:
        print("\nThen: python fetch_inputs.py --inputs-dir "
              f"{args.inputs_dir} --verify")


if __name__ == "__main__":
    main()
