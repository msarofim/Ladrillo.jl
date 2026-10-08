"""
forcing_inputs.py -- the ONE table of every input file the forcing builders read.

Each logical input has
  * `rel`       its path under an --inputs-dir, as fetch_inputs.py lays it out
                (and as inputs.sha256 lists it);
  * `fairtofredi`  where the research repo keeps the same file, so a build can be
                pointed at the copies already in the maintainer's research
                checkout with --fairtofredi-root instead of a fetched --inputs-dir;
  * `public`    True / False, and `source` / `url` / `licence` saying from where.
                `public=False` means a reviewer CANNOT fetch it: it is gated,
                private, or of unestablished provenance. Nothing is hidden behind
                a tolerance here either. A build without such a file skips the
                steps that read it (build_forcing.py says which files it then
                cannot build); the README says which ones nobody outside can get.

The three RCMIP files are special: fair's `fill_from_rcmip()` reads them through
pooch's cache, not through a path we pass, so they always resolve to that cache
(see `rcmip_cache_path`). fetch_inputs.py primes the cache and also keeps a copy
under <inputs-dir>/rcmip/ so inputs.sha256 can be checked with `shasum -c`.

Nothing in this module downloads anything.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHA_FILE = HERE / "inputs.sha256"

# -- Named constants ----------------------------------------------------------
ZENODO_CALIB_RECORD = "18828694"     # fair-calibrate v1.6.0 calibration data
ZENODO_CALIB_LABEL = ("Zenodo 10.5281/zenodo.18828694 'fair calibration data' "
                      "v1.6.0 (2026-03-02)")
ZENODO_CALIB_URL = "https://zenodo.org/records/{rec}/files/{name}"
ZENODO_VV_RECORD = "20713982"        # van Vuuren ScenarioMIP extensions, v1.1.1
ZENODO_VV_LABEL = ("Zenodo 10.5281/zenodo.20713982 v1.1.1 (2026-06-16), "
                   "chrisroadmap/cmip7-scenariomip-1.1.1.zip")
ZENODO_VV_ZIP_URL = ("https://zenodo.org/api/records/20713982/files/"
                     "chrisroadmap/cmip7-scenariomip-1.1.1.zip/content")
ZENODO_VV_MEMBER_SUFFIX = "data/fair-inputs/emissions_1750-2500.csv"

# Licences, read from the Zenodo and GitHub APIs (metadata only) on 2026-10-08.
ZENODO_CALIB_LICENCE = "CC-BY-4.0 (Zenodo API, 2026-10-08); attribution required"
ZENODO_CALIB_ATTRIBUTION = ("Smith, C. (2026). fair calibration data (v1.6.0) "
                            "[Data set]. Zenodo. https://doi.org/10.5281/zenodo.18828694")
ZENODO_VV_LICENCE = "MIT (Zenodo API, 2026-10-08)"

# chrisroadmap/fair-calibrate, pinned to the last commit that touched the file,
# as reported by the GitHub API (commits?path=..., 2026-10-08; metadata only,
# nothing downloaded). Whether that commit's bytes are the copy the shipped
# forcing was built from is NOT yet verified: fetch_inputs.py --verify decides
# on first use.
FAIRCAL_CHRIS_COMMIT = "3f72545665f61c0416299412be68e4deab46ec61"
FAIRCAL_CHRIS_PATH = "data/emissions/CH4Pledge-rcmip-emissions-ssp245-CLE-MFR-filled.csv"
FAIRCAL_CHRIS_URL = ("https://raw.githubusercontent.com/chrisroadmap/fair-calibrate/"
                     f"{FAIRCAL_CHRIS_COMMIT}/{FAIRCAL_CHRIS_PATH}")
FAIRCAL_LICENCE = "MIT (GitHub API, chrisroadmap/fair-calibrate, 2026-10-08)"

# fair 2.2.4's own RCMIP v5.1.0 retrievals (fair/io/fill_from.py), copied so the
# offline guard checks EXACTLY what fill_from_rcmip() will ask pooch for.
RCMIP_BASE = "https://rcmip-protocols-au.s3-ap-southeast-2.amazonaws.com/v5.1.0/"
RCMIP_FILES = {
    "rcmip_emissions": ("rcmip-emissions-annual-means-v5-1-0.csv",
                        "md5:4044106f55ca65b094670e7577eaf9b3"),
    "rcmip_concentrations": ("rcmip-concentrations-annual-means-v5-1-0.csv",
                             "md5:0d82c3c3cdd4dd632b2bb9449a5c315f"),
    "rcmip_forcing": ("rcmip-radiative-forcing-annual-means-v5-1-0.csv",
                      "md5:87ef6cd4e12ae0b331f516ea7f82ccba"),
}


def _z(name, what):
    return dict(rel=f"zenodo_{ZENODO_CALIB_RECORD}/{name}",
                fairtofredi=f"v160_cache/{name}", public=True,
                url=ZENODO_CALIB_URL.format(rec=ZENODO_CALIB_RECORD, name=name),
                source=f"{ZENODO_CALIB_LABEL}, file {name}",
                licence=ZENODO_CALIB_LICENCE, attribution=ZENODO_CALIB_ATTRIBUTION,
                what=what)


INPUTS = {
    # -- fair-calibrate 1.6.0 (Zenodo 18828694) ------------------------------
    "calib_params": _z("calibrated_constrained_parameters.csv",
                       "841-config posterior"),
    "calib_species": _z("species_configs_properties.csv",
                        "62-species configs (stock 1.6.0)"),
    "calib_volcanic": _z("volcanic_forcing_timebounds_cmip7.csv", "CMIP7 volcanic ERF"),
    "calib_solar": _z("solar_forcing_timebounds_cmip7.csv", "CMIP7 solar ERF"),
    "calib_irrigation": _z("irrigation_forcing_timebounds_cmip7.csv",
                           "per-marker irrigation ERF"),
    "calib_landuse": _z("land_use_forcing_timebounds_cmip7.csv",
                        "per-marker land-use ERF"),
    "cmip7_hist_emissions": _z("historical_emissions_1750-2023_cmip7.csv",
                               "CMIP7 historical emissions 1750-2023"),
    # -- RCMIP v5.1.0 (read by fair.fill_from_rcmip through pooch) -----------
    **{k: dict(rel=f"rcmip/{fn}", fairtofredi=None, public=True,
               url=RCMIP_BASE + fn, md5=md5,
               source="RCMIP v5.1.0 (Nicholls et al. 2020), fetched by fair 2.2.4 "
                      "fill_from_rcmip() via pooch",
               licence="RCMIP terms (NOT checked)",
               what="RCMIP SSP future emissions (only the emissions file's values "
                    "are used; fair retrieves all three)")
       for k, (fn, md5) in RCMIP_FILES.items()},
    # -- SSP2-4.5 (calibration forcing) future and its post-2100 shape -------
    "chris_ssp245": dict(
        rel="chrisroadmap_fair-calibrate/CH4Pledge-rcmip-emissions-ssp245-CLE-MFR-filled.csv",
        fairtofredi="calibration_v145/raw_chrisroadmap_ssp245/ssp245_chrisroadmap.csv",
        public=True, url=FAIRCAL_CHRIS_URL,
        source=f"GitHub chrisroadmap/fair-calibrate @ {FAIRCAL_CHRIS_COMMIT[:7]}, "
               f"{FAIRCAL_CHRIS_PATH} (MESSAGE-GLOBIOM SSP2-4.5 + CH4-pledge "
               "variants). Commit from GitHub API metadata; the byte identity of "
               "the fetched file is UNVERIFIED (no download was made) -- the "
               "sha256 check decides on first use.",
        licence=FAIRCAL_LICENCE,
        what="SSP2-4.5 future 2024-2100"),
    "ar6_ssp245": dict(
        rel="edf_rennels/AR6_emissions_ssp245_1750_2300.csv",
        fairtofredi="AR6_emissions_ssp245_1750_2300.csv",
        public=False, url=None,
        source="EDF / Lisa Rennels' AR6 SSP2-4.5 1750-2300 file, as received in the "
               "research repo (2026-06-01). No public origin is recorded for it.",
        licence="NONE RECORDED",
        what="post-2100 r(t) SHAPE for the SSP2-4.5 extension (ratios only)"),
    # -- van Vuuren ScenarioMIP-CMIP7 markers --------------------------------
    "vv_ext_zenodo": dict(
        rel=f"zenodo_{ZENODO_VV_RECORD}/emissions_1750-2500.csv",
        fairtofredi="data/vanvuuren/zenodo_v111_emissions_1750-2500.csv",
        public=True, url=ZENODO_VV_ZIP_URL, zip_member_suffix=ZENODO_VV_MEMBER_SUFFIX,
        source=f"{ZENODO_VV_LABEL}, member {ZENODO_VV_MEMBER_SUFFIX}. NOT the DOI "
               "printed in van Vuuren et al. 2026 (10.5281/zenodo.14382495 = "
               "stylized v1.0).",
        licence=ZENODO_VV_LICENCE,
        what="every van Vuuren value from 2023.5 on (the 2023.5 anchor and the "
             "2024.5-2300.5 scenario)"),
}


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def md5(path, chunk=1 << 20):
    h = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def rcmip_cache_path(key):
    """Where pooch.retrieve(url, known_hash) -- fair's call -- keeps the file."""
    import pooch
    from pooch.utils import unique_file_name
    fn, _ = RCMIP_FILES[key]
    return Path(pooch.os_cache("pooch")) / unique_file_name(RCMIP_BASE + fn)


def require_rcmip_cache():
    """Exit unless every RCMIP file fair will retrieve is cached with fair's md5.

    pooch.retrieve() downloads whatever is missing or fails its hash; this check
    runs first so the build never touches the network.
    """
    bad = []
    for key, (fn, want) in RCMIP_FILES.items():
        p = rcmip_cache_path(key)
        if not p.exists():
            bad.append(f"{fn}: not in pooch cache ({p})")
        elif "md5:" + md5(p) != want:
            bad.append(f"{fn}: cached copy fails fair's {want}")
    if bad:
        sys.exit("ERROR: RCMIP v5.1.0 is not fully cached, so fill_from_rcmip() "
                 "would DOWNLOAD it:\n  " + "\n  ".join(bad)
                 + "\nRun fetch_inputs.py (it primes pooch's cache), then retry.")
    print("  [RCMIP-CACHE] all three RCMIP v5.1.0 files cached with fair's md5  PASS")


def expected_sha256():
    """{rel: sha256} from inputs.sha256."""
    out = {}
    for line in SHA_FILE.read_text().splitlines():
        if line.strip():
            h, rel = line.split(None, 1)
            out[rel.strip()] = h
    return out


def resolve(key, inputs_dir=None, fairtofredi_root=None):
    """Absolute path of one logical input, from exactly one of the two roots."""
    spec = INPUTS[key]
    if key in RCMIP_FILES:
        return rcmip_cache_path(key)
    if (inputs_dir is None) == (fairtofredi_root is None):
        sys.exit("ERROR: give exactly one of --inputs-dir / --fairtofredi-root")
    if inputs_dir is not None:
        return (Path(inputs_dir) / spec["rel"]).resolve()
    return (Path(fairtofredi_root).expanduser() / spec["fairtofredi"]).resolve()


def check_inputs(inputs_dir=None, fairtofredi_root=None, keys=None, quiet=False):
    """{key: (status, detail)} with status "ok", "absent" or "bad".

    "absent" = the file is not there (for RCMIP: not in pooch's cache). "bad" =
    it is there and fails inputs.sha256 (RCMIP: or fair's own md5), or the table
    does not list it. build_forcing.py skips what needs an absent input and
    ABORTS on any bad one: a wrong file is never built on.
    """
    want = expected_sha256()
    out = {}
    for key in keys or INPUTS:
        rel = INPUTS[key]["rel"]
        p = resolve(key, inputs_dir, fairtofredi_root)
        if rel not in want:
            out[key] = ("bad", "not listed in inputs.sha256")
        elif not p.exists():
            out[key] = ("absent", f"missing: {p}")
        elif sha256(p) != want[rel]:
            out[key] = ("bad", f"sha256 {sha256(p)} != {want[rel]}  ({p})")
        elif key in RCMIP_FILES and "md5:" + md5(p) != RCMIP_FILES[key][1]:
            out[key] = ("bad", f"cached copy fails fair's {RCMIP_FILES[key][1]}  ({p})")
        else:
            out[key] = ("ok", str(p))
        if not quiet:
            st, why = out[key]
            print(f"  {st.upper():<6} {key:<22} {rel}"
                  + ("" if st == "ok" else f"\n         {why}"))
    return out


def verify(inputs_dir=None, fairtofredi_root=None, keys=None, quiet=False):
    """Hash every resolved input against inputs.sha256. Returns the failures."""
    want = expected_sha256()
    fails = []
    for key in keys or INPUTS:
        rel = INPUTS[key]["rel"]
        p = resolve(key, inputs_dir, fairtofredi_root)
        if rel not in want:
            fails.append((key, "not listed in inputs.sha256"))
            continue
        if not p.exists():
            fails.append((key, f"missing: {p}"))
            continue
        got = sha256(p)
        ok = got == want[rel]
        if not ok:
            fails.append((key, f"sha256 {got} != {want[rel]}  ({p})"))
        if not quiet:
            print(f"  {'PASS' if ok else 'FAIL'}  {key:<22} {rel}")
    return fails
