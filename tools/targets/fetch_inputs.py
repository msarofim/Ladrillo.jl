#!/usr/bin/env python3
"""
fetch_inputs.py -- lay out the inputs Ladrillo does NOT ship in one
--inputs-dir, and check them against inputs.sha256.

    python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs            # the plan, no network
    python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --download # Frederikse (2020), Zenodo
    python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --verify

Layout of --inputs-dir: each `fetch/<rel>` of inputs.sha256 at <inputs-dir>/<rel>,
each `raw/<rel>` at <inputs-dir>/raw/<rel>.

The standard rebuild needs only the two Frederikse (2020) files: everything else
is shipped in data/sources/. --download fetches them from Zenodo 10621070 (their
md5 there equals the copies the targets were built from).

CURRENT VINTAGES (optional). The NOAA STAR and NCEI files and the Berkeley Earth
and GISTEMP Greenland series in data/sources/ come from rolling products. To
build on today's versions instead: --download also writes the current NOAA STAR
file, reduced to annual means exactly as the research did (download_obs.py), to
<inputs-dir>/noaa/nasa_gmsl_annual.csv; place a current NCEI file at
<inputs-dir>/noaa/noaa_thermosteric_w0-2000m_yearly.dat; put current Berkeley
Earth and GISTEMP grids (and the GTN-G shapefile) in <inputs-dir>/raw/ and run
--extract --no-pin-check. Then build with build_targets.py --current-vintages.
Today's vintages do NOT reproduce the shipped bytes.

--from-research ROOT (maintainer) copies the fetch files out of an SLR-RFF-BRICK
checkout and symlinks its raw products into <inputs-dir>/raw/.

--extract runs extract/extract_greenland_zones.py --products BerkeleyEarth
GISTEMP on <inputs-dir>/raw/ and writes <inputs-dir>/greenland_zones/. Without
--no-pin-check the raw files must match their pins (the archived vintages), and
the output then equals the shipped files.
"""
from __future__ import annotations

import argparse
import io
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import target_inputs as TI  # noqa: E402

UA = "Mozilla/5.0 (Ladrillo tools/targets fetch_inputs.py)"
STAR_URLS = ("https://www.star.nesdis.noaa.gov/socd/lsa/SeaLevelRise/slr/slr_sla_gbl_free_ref_90.csv",
             "https://www.star.nesdis.noaa.gov/socd/lsa/SeaLevelRise/slr/slr_sla_gbl_keep_ref_90.csv")
STAR_REL = "noaa/nasa_gmsl_annual.csv"
FREDERIKSE = {"frederikse2020/frederikse2020_global_basin_timeseries.xlsx": "global_basin_timeseries.xlsx",
              "frederikse2020/frederikse2020_GMSL_ensembles.nc": "GMSL_ensembles_F20.nc"}
GREENLAND_FETCH_PRODUCTS = ["BerkeleyEarth", "GISTEMP"]
CURRENT_VINTAGE_KEYS = {"noaa_star", "noaa_steric", "greenland_berkeley", "greenland_gistemp"}


def http_get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def star_annual(text: str):
    """NOAA STAR GMSL csv -> annual means: research download_obs.py fetch_nasa_gmsl(), verbatim."""
    import numpy as np
    import pandas as pd
    decimal_years, slas = [], []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("HDR") or s.startswith("#") or s.startswith("%"):
            continue
        parts = [p.strip() for p in s.replace(",", " ").split()]
        if len(parts) < 2:
            continue
        try:
            dy = float(parts[0])
            sla = float(parts[1])
        except ValueError:
            continue
        if dy < 1990 or dy > 2050:
            continue
        decimal_years.append(dy)
        slas.append(sla)
    df = pd.DataFrame({"decimal_year": decimal_years, "sla_mm": slas})
    df["year"] = df["decimal_year"].astype(int)
    annual = (df.groupby("year").agg(value=("sla_mm", "mean"), n=("sla_mm", "size")).reset_index())
    annual = annual[annual["n"] >= 20].drop(columns="n")
    annual["sigma"] = np.nan
    annual = annual.sort_values("year").reset_index(drop=True)
    buf = io.StringIO()
    annual[["year", "value", "sigma"]].to_csv(buf, index=False)
    return buf.getvalue().encode()


def status(inputs: Path, pins: dict):
    rows = []
    for k, sha in pins.items():
        kind, rel = k.split("/", 1)
        if kind == "sources":
            continue
        p = inputs / rel if kind == "fetch" else inputs / "raw" / rel
        st = "absent" if not p.exists() else ("OK" if TI.sha256(p) == sha else "MISMATCH")
        rows.append((kind, rel, st))
    return rows


def plan(inputs: Path, pins: dict):
    st = {(k, r): s for k, r, s in status(inputs, pins)}
    print(f"inputs-dir {inputs}\n")
    print("FETCH YOURSELF (not shipped: redistribution in doubt; ruling 2026-11-12). Needed by every build;"
          " --download gets them.")
    for s in TI.SOURCES:
        for kind, rel, _ in s["files"]:
            if kind != "fetch":
                continue
            print(f"  [{st.get(('fetch', rel), 'absent'):8s}] {rel}")
            print(f"      {s['short']}: {s['name']}")
            print(f"      from: {s['url']}")
            print(f"      licence: {s['licence']}" + ("   ROLLING: a fresh download will not match" if s["rolling"] else ""))
    print("\nCURRENT VINTAGES (optional; not needed). Shipped files from rolling products; today's versions do NOT")
    print("reproduce the shipped bytes. Build on them with build_targets.py --current-vintages.")
    for s in TI.SOURCES:
        for kind, rel, _ in s["files"]:
            if kind == "sources" and s["rolling"] and s["key"] in CURRENT_VINTAGE_KEYS:
                p = inputs / rel
                print(f"  [{'present' if p.exists() else 'absent':8s}] {rel}   {s['short']}: {s['url']}")
    print("\nRAW (maintainer re-extraction only; --extract needs the Berkeley, GISTEMP and GTN-G files).")
    for r in TI.RAW:
        print(f"  [{st.get(('raw', r['rel']), 'absent'):8s}] raw/{r['rel']}   {r['url']}"
              + (f"   ({r['note']})" if r["note"] else ""))


def from_research(root: Path, inputs: Path):
    for kind, rel, research, _ in TI.files("fetch"):
        if research is None:
            continue
        dst = inputs / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / research, dst)
        print(f"  copied  {rel}")
    for r in TI.RAW:
        dst = inputs / "raw" / r["rel"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(root / TI.RESEARCH_RAW / r["rel"])
        print(f"  linked  raw/{r['rel']}")


def download(inputs: Path):
    for rel, name in FREDERIKSE.items():
        url = TI.ZENODO_DANG.format(name=name)
        print(f"  GET {url}")
        dst = inputs / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(http_get(url))
    for url in STAR_URLS:
        print(f"  GET {url}  (ROLLING: the pinned annual means are from the 2026-05-20 release)")
        try:
            raw = http_get(url, timeout=60)
        except OSError as e:
            print(f"    failed: {e}")
            continue
        if len(raw) > 500:
            dst = inputs / STAR_REL
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(star_annual(raw.decode("utf-8", errors="replace")))
            break
    print(f"  {STAR_REL} is the CURRENT STAR release: it does not reproduce the shipped 2026-05-20 file and is used "
          "only by build_targets.py --current-vintages")
    print("  not downloaded (see the plan): NOAA NCEI thermosteric, and every raw/ product")


def extract(inputs: Path, pin_check=True):
    cmd = [sys.executable, str(HERE / "extract" / "extract_greenland_zones.py"), "--raw-dir", str(inputs / "raw"),
           "--out-dir", str(inputs / "greenland_zones"), "--products", *GREENLAND_FETCH_PRODUCTS]
    if not pin_check:
        cmd.append("--no-verify")
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    subprocess.run(cmd, check=True, env=env)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs-dir", required=True)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--from-research", metavar="ROOT")
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--no-pin-check", action="store_true", help="with --extract: allow current (unpinned) raw grids")
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args()
    inputs = Path(a.inputs_dir).expanduser().resolve()
    inputs.mkdir(parents=True, exist_ok=True)
    pins = TI.read_pins()
    if a.from_research:
        from_research(Path(a.from_research).expanduser().resolve(), inputs)
    if a.download:
        download(inputs)
    if a.extract:
        extract(inputs, pin_check=not a.no_pin_check)
    if a.verify:
        bad = 0
        for kind, rel, st in status(inputs, pins):
            if kind == "fetch" or st != "absent":
                print(f"  {st:8s} {kind}/{rel}")
                bad += st != "OK" and kind == "fetch"
        print("VERIFY: " + ("all fetch inputs present and pinned" if not bad else f"{bad} fetch input(s) absent or mismatched"))
        sys.exit(1 if bad else 0)
    if not (a.from_research or a.download or a.extract):
        plan(inputs, pins)


if __name__ == "__main__":
    main()
