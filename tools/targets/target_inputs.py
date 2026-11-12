"""
target_inputs.py -- the ONE table of every input the calibration-target and
driver builders read, where it comes from, and under what licence.

Kinds (the prefix each file carries in inputs.sha256):

  sources/  SHIPPED in data/sources/. The licence clearly allows
            redistribution (attribution at most). Copied byte-for-byte from the
            research repository (SLR-RFF-BRICK), or a small extraction of a
            large product (its builder is in tools/targets/extract/).
  fetch/    FETCH YOURSELF into an --inputs-dir, at the path shown. Not shipped
            because the licence leaves redistribution in doubt (Marcus's ruling
            of 2026-11-12: such files are not shipped; they are pinned by sha256
            and the builder refuses any other file). Some are ROLLING products
            the publisher overwrites: a download made today will not match the
            pin (see `rolling`).
  raw/      RAW gridded products read only by the maintainer re-extractions
            (build_targets.py --from-raw, fetch_inputs.py --extract).

Licences were read from the publishers' pages and the Zenodo API on
2026-11-12 (metadata only; nothing was downloaded).

Nothing in this module downloads anything.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHA_FILE = HERE / "inputs.sha256"
REPO_ROOT = HERE.parent.parent
SOURCES_DIR = REPO_ROOT / "data" / "sources"
LICENCES_READ = "2026-11-12"
ZENODO_DANG = "https://zenodo.org/records/10621070/files/{name}?download=1"

CMIP6_MODELS = [  # the 40 files of research data/cmip6_gis/, in its sorted order
    "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CAMS-CSM1-0", "CAS-ESM2-0", "CESM2-WACCM",
    "CESM2", "CIESM", "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1-HR", "CNRM-CM6-1", "CNRM-ESM2-1",
    "CanESM5-CanOE", "CanESM5", "EC-Earth3-Veg-LR", "EC-Earth3-Veg", "EC-Earth3", "FGOALS-f3-L", "FGOALS-g3",
    "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "IITM-ESM", "INM-CM4-8",
    "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "KIOST-ESM", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR",
    "MPI-ESM1-2-LR", "MRI-ESM2-0", "NESM3", "NorESM2-LM", "NorESM2-MM",
]

GLAMBIE_FILES = ["0_global.csv", "1_alaska.csv", "2_western_canada_us.csv", "3_arctic_canada_north.csv",
                 "4_arctic_canada_south.csv", "5_greenland_periphery.csv", "6_iceland.csv", "7_svalbard.csv",
                 "8_scandinavia.csv", "9_russian_arctic.csv", "10_north_asia.csv", "11_central_europe.csv",
                 "12_caucasus_middle_east.csv", "13_central_asia.csv", "14_south_asia_west.csv",
                 "15_south_asia_east.csv", "16_low_latitudes.csv", "17_southern_andes.csv", "18_new_zealand.csv",
                 "19_antarctic_and_subantarctic.csv"]

# One entry per source dataset. `files`: (kind, rel, research path or None).
SOURCES = [
    dict(key="frederikse2020", short="Frederikse et al. (2020)",
         name="sea-level budget: global basin time series (xlsx) and the 5000-member weighted ensemble (nc)",
         files=[("fetch", "frederikse2020/frederikse2020_global_basin_timeseries.xlsx",
                 "data/observations/raw/frederikse2020_global_basin_timeseries.xlsx"),
                ("fetch", "frederikse2020/frederikse2020_GMSL_ensembles.nc",
                 "data/observations/raw/frederikse2020_GMSL_ensembles.nc")],
         url=ZENODO_DANG.format(name="global_basin_timeseries.xlsx") + " and "
             + ZENODO_DANG.format(name="GMSL_ensembles_F20.nc") + " (rename on save)",
         doi="10.5281/zenodo.10621070 (redistributed by Dangendorf et al.; original 10.5281/zenodo.3862995)",
         accessed="2026-05-21 (xlsx), 2026-07-19 (nc)", rolling=False,
         licence="original record 3862995: 'Other (Not Open)', all rights reserved; record 10621070: CC BY 4.0, "
                 "which cannot relicense Frederikse's files",
         verified="local copies match Zenodo 10621070's published md5 (API, 2026-11-12)",
         used_by="every component target (ais, gis, gsic, steric, lws) 1900-2018 and its sigma; the total's "
                 "budget-closure sigma; the land-water splice offset"),
    dict(key="dangendorf2024_v2", short="Dangendorf et al. (2024), corrected global file",
         name="Kalman-smoother GMSL reconstruction, global series with standard errors (v2)",
         files=[("sources", "dangendorf2024/dangendorf2024_KalmanSmootherHR_Global_v2.nc",
                 "data/observations/raw/dangendorf2024_KalmanSmootherHR_Global_v2.nc")],
         url="", doi="none: personal communication, S. Dangendorf to M. Sarofim, 2026-08-07 (fixes the slot "
                     "shift in Zenodo 10621070's KalmanSmootherHR_Global.nc)",
         accessed="2026-08-07", rolling=False,
         licence="none stated; shipped by Marcus's ruling (2026-11-12)", verified="",
         used_by="dang_sig, dang_lo, dang_hi 1900-2021 (the total target's own SE)"),
    dict(key="dangendorf2024_gmsl", short="Dangendorf et al. (2024)",
         name="annual GMSL = cos(lat)-weighted mean of the HR field (extraction)",
         files=[("sources", "dangendorf2024/dangendorf2024_gmsl_annual.csv",
                 "data/observations/dangendorf2024_gmsl_annual.csv")],
         url="https://zenodo.org/records/10621070", doi="10.5281/zenodo.10621070",
         accessed="2026-07-20", rolling=False, licence="CC BY 4.0 (Zenodo API)",
         verified="the two raw files it is extracted from match Zenodo's md5", extractor="extract/extract_dangendorf_gmsl.py",
         used_by="dang (the total target) 1900-2021"),
    dict(key="grace_mass", short="JPL GRACE/GRACE-FO mascons RL06.3Mv04",
         name="Antarctic and Greenland mass time series",
         files=[("sources", "grace/grace_antarctica_mass.txt", "data/observations/raw/grace_antarctica_mass.txt"),
                ("sources", "grace/grace_greenland_mass.txt", "data/observations/raw/grace_greenland_mass.txt")],
         url="https://podaac.jpl.nasa.gov/ (Earthdata login)", doi="10.5067/TEMSC-AT634, 10.5067/TEMSC-GT634",
         accessed="2026-06-13", rolling=True,
         licence="NASA EOSDIS data-use guidance: openly shared without restriction (CC0 unless marked)",
         verified="", used_by="ais and gis 2019-2026; the mask gates of the land-water extraction"),
    dict(key="grace_lws", short="JPL GRACE/GRACE-FO mascons RL06.3Mv04",
         name="monthly mass over global land minus both ice sheets (extraction)",
         files=[("sources", "grace/grace_lws_domain_monthly.csv", None)],
         url="https://podaac.jpl.nasa.gov/ (Earthdata login)", doi="10.5067/TEMSC-3JC634",
         accessed="2026-08-24", rolling=True, licence="as above (NASA, without restriction)", verified="",
         extractor="extract/extract_grace_lws_domain.py", used_by="lws 2019-2023, held flat to 2026"),
    dict(key="glambie", short="GlaMBIE Team (2025)",
         name="regional calendar-year glacier mass change (results, calendar_years/ of glambie_data.zip)",
         files=[("sources", f"glambie/calendar_years/{n}",
                 f"glambie_data.zip:glambie_results_20240716/calendar_years/{n}") for n in GLAMBIE_FILES],
         url="https://doi.org/10.5904/wgms-glambie-2024-07", doi="10.5904/wgms-glambie-2024-07",
         accessed="2026-06-13", rolling=False, licence="WGMS data policy: open access with citation (= CC BY 4.0)",
         verified="", used_by="gsic 2019-2023 (global minus region 5); glacier block weights, modern rates, 2020 "
                             "losses; the region-19 seam; glaciers removed from the land-water extension"),
    dict(key="noaa_steric", short="NOAA NCEI",
         name="World Ocean 0-2000 m thermosteric sea level, yearly",
         files=[("fetch", "noaa/noaa_thermosteric_w0-2000m_yearly.dat",
                 "data/observations/raw/noaa_thermosteric_w0-2000m_yearly.dat")],
         url="https://www.ncei.noaa.gov/access/global-ocean-heat-content/ (yearly thermosteric sea level, "
             "World Ocean, 0-2000 m)", doi="", accessed="2026-06-13", rolling=True,
         licence="none stated on the pages read", verified="",
         used_by="steric 2019-2025"),
    dict(key="noaa_star", short="NOAA STAR Laboratory for Satellite Altimetry",
         name="global mean sea level, seasonal signals removed; annual means of the 2026-05 release",
         files=[("fetch", "noaa/nasa_gmsl_annual.csv", "data/observations/nasa_gmsl_annual.csv")],
         url="https://www.star.nesdis.noaa.gov/socd/lsa/SeaLevelRise/slr/slr_sla_gbl_free_ref_90.csv "
             "(fetch_inputs.py --download derives the annual means as the research did)",
         doi="", accessed="2026-05-20", rolling=True,
         licence="no licence; LSA data-use policy (2008) requires a credit line", verified="",
         used_by="dang 2022-2024"),
    dict(key="gmip3_tables", short="GlacierMIP3 (Zekollari et al. 2025)",
         name="Tables S1a and S3, region characteristics",
         files=[("sources", "gmip3/table_S1a.csv", "data/observations/raw/gmip3/table_S1a.csv"),
                ("sources", "gmip3/table_S3.csv", "data/observations/raw/gmip3/table_S3.csv"),
                ("sources", "gmip3/3_shift_summary_region_characteristicsFeb12_2024.csv",
                 "data/observations/raw/gmip3/3_shift_summary_region_characteristicsFeb12_2024.csv")],
         url="https://zenodo.org/records/15046588", doi="10.5281/zenodo.15046588", accessed="2026-08-07",
         rolling=False, licence="CC BY 4.0 (Zenodo API)", verified="",
         used_by="glacier inventory, response times, ISIMIP3 amplification"),
    dict(key="gmip3_ladder", short="GlacierMIP3 (Zekollari et al. 2025)",
         name="committed-loss ladder per Ladrillo glacier reservoir (extraction)",
         files=[("sources", "gmip3/gmip3_block_ladder_L27.csv", "outputs/d1_block_ladder_cache.csv")],
         url="https://zenodo.org/records/15046588", doi="10.5281/zenodo.15046588", accessed="2026-08-07",
         rolling=False, licence="CC BY 4.0 (Zenodo API)", verified="",
         extractor="research python/d1_multireservoir_cell.py block_ladder() (not ported; see README)",
         used_by="the glacier equilibrium curves (b, T_off) and their committed rungs"),
    dict(key="hadcrut5_glacier", short="Met Office HadCRUT5.0.2.0 (Morice et al. 2021)",
         name="glacier-area-weighted annual temperature per RGI region, and the global mean (extraction)",
         files=[("sources", "hadcrut5/t_glac_hadcrut5.csv", "data/observations/t_glac_hadcrut5.csv"),
                ("sources", "hadcrut5/t_glac_regions_hadcrut5.csv", "data/observations/t_glac_regions_hadcrut5.csv")],
         url="https://www.metoffice.gov.uk/hadobs/hadcrut5/", doi="10.1029/2019JD032361", accessed="2026-08-06",
         rolling=True, licence="Open Government Licence v3 (Crown copyright); acknowledgement line required",
         verified="", extractor="extract/build_t_glac.py",
         used_by="the three glacier reservoir drivers; the observational amplification fit"),
    dict(key="greenland_hadcrut5", short="Met Office HadCRUT5.0.2.0",
         name="Greenland land-zone and global annual temperature (extraction)",
         files=[("sources", "greenland_zones/greenland_zone_series_HadCRUT5.csv", None)],
         url="https://www.metoffice.gov.uk/hadobs/hadcrut5/", doi="10.1029/2019JD032361", accessed="2026-08-06",
         rolling=True, licence="Open Government Licence v3 (the zone weights use the Berkeley Earth 1-degree land "
                               "fraction and the GTN-G region-05 polygon; no Berkeley values are in the file)",
         verified="", extractor="extract/extract_greenland_zones.py",
         used_by="the Greenland driver; the amplification prior and anchor"),
    dict(key="greenland_berkeley", short="Berkeley Earth",
         name="Greenland land-zone and global annual temperature (extraction, made by fetch_inputs.py --extract)",
         files=[("fetch", "greenland_zones/greenland_zone_series_BerkeleyEarth.csv", None)],
         url="https://berkeleyearth.org/data/ (Land_and_Ocean_LatLong1.nc)", doi="", accessed="2026-08-09",
         rolling=True, licence="CC BY-NC 4.0 (berkeleyearth.org/data)", verified="",
         extractor="extract/extract_greenland_zones.py --products BerkeleyEarth",
         used_by="the amplification prior and anchor (one of three products)"),
    dict(key="greenland_gistemp", short="NASA GISS GISTEMP v4",
         name="Greenland land-zone and global annual temperature (extraction, made by fetch_inputs.py --extract)",
         files=[("fetch", "greenland_zones/greenland_zone_series_GISTEMP.csv", None)],
         url="https://data.giss.nasa.gov/gistemp/ (gistemp1200_GHCNv4_ERSSTv5.nc.gz)", doi="", accessed="2026-08-09",
         rolling=True, licence="none stated (citation requested)", verified="",
         extractor="extract/extract_greenland_zones.py --products GISTEMP",
         used_by="the amplification prior and anchor (one of three products)"),
    dict(key="cmip6_gis", short="CMIP6 (40 models)",
         name="Amon tas, annual global and Greenland-zone means, historical + ssp126/245/585 (extraction)",
         files=[("sources", f"cmip6_gis/tas_series_gis_{m}.csv", f"data/cmip6_gis/tas_series_gis_{m}.csv")
                for m in CMIP6_MODELS],
         url="https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv (Pangeo)", doi="",
         accessed="2026-08-12 / 2026-08-18", rolling=True,
         licence="CC BY 4.0 for all groups since 2022-10 (CMIP6 Terms of Use v6.2); GISS sources CC0",
         verified="", extractor="extract/reduce_cmip6_tas_gis.py (provenance; not re-run)",
         used_by="the Greenland amplification shape S(dT)"),
]

# Raw products: read only by the maintainer re-extractions.
RAW = [
    dict(rel="HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc", accessed="2026-08-06", mb=32.5,
         url="https://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.0.2.0/analysis/"
             "HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc",
         licence="OGL v3", note="updated monthly in place"),
    dict(rel="Land_and_Ocean_LatLong1.nc", accessed="2026-08-09", mb=454.6,
         url="https://berkeleyearth.org/data/ (Gridded, Land + Ocean, 1x1 degree)",
         licence="CC BY-NC 4.0", note="updated in place"),
    dict(rel="gistemp1200_GHCNv4_ERSSTv5.nc", accessed="2026-08-09", mb=57.0,
         url="https://data.giss.nasa.gov/pub/gistemp/gistemp1200_GHCNv4_ERSSTv5.nc.gz (gunzip)",
         licence="none stated", note="updated monthly in place"),
    dict(rel="GRCTellus.JPL.200204_202606.GLO.RL06.3M.MSCNv04CRI.nc", accessed="2026-08-24", mb=45.8,
         url="https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4 (Earthdata login)",
         licence="NASA, without restriction", note="the name carries the last month; replaced monthly"),
    dict(rel="GlacReg_2023/GTN-G_202307_o1regions.shp", accessed="2026-08-06", mb=None,
         url="https://doi.org/10.5904/gtng-glacreg-2023-07 (GlacReg_2023.zip)",
         licence="GTN-G data policy: open access with citation (= CC BY 4.0), also 'public domain for "
                 "non-commercial use'", note=""),
    dict(rel="GlacReg_2023/GTN-G_202307_o1regions.shx", accessed="2026-08-06", mb=None, url="(same zip)",
         licence="as above", note=""),
    dict(rel="GlacReg_2023/GTN-G_202307_o1regions.dbf", accessed="2026-08-06", mb=None, url="(same zip)",
         licence="as above", note=""),
    dict(rel="dangendorf2024_KalmanSmootherHR_Fields.nc", accessed="2026-07-20", mb=729.5,
         url=ZENODO_DANG.format(name="KalmanSmootherHR_Fields.nc"), licence="CC BY 4.0",
         note="matches Zenodo's md5"),
    dict(rel="dangendorf2024_KalmanSmootherHR_Global.nc", accessed="2026-07-19", mb=0.06,
         url=ZENODO_DANG.format(name="KalmanSmootherHR_Global.nc"), licence="CC BY 4.0",
         note="matches Zenodo's md5; slot-shifted upstream, read only for the indicative-SE column"),
    dict(rel="dangendorf2024_grid_lonlat.csv", accessed="2026-07-20", mb=1.2,
         url="(lon/lat of the Fields file's grid points; how it was produced is not recorded in the research repo)",
         licence="CC BY 4.0 (derived from Zenodo 10621070)", note="origin not recorded"),
]
RESEARCH_RAW = "data/observations/raw"      # where the research checkout keeps every RAW file


def files(kind=None):
    """[(kind, rel, research path or None, source key)]."""
    return [(k, rel, res, s["key"]) for s in SOURCES for (k, rel, res) in s["files"] if kind in (None, k)]


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_pins(path=SHA_FILE):
    """{'<kind>/<rel>': sha256} from inputs.sha256."""
    pins = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        sha, rel = line.split(None, 1)
        pins[rel.strip()] = sha
    return pins


KIND_LABEL = {"sources": "shipped", "fetch": "fetch yourself"}


def render_markdown() -> str:
    """The README sources table, rendered from SOURCES and RAW (so the two cannot drift)."""
    out = ["| source | what | files | kind | where from | accessed | licence (read " + LICENCES_READ + ") | feeds |",
           "|---|---|---|---|---|---|---|---|"]
    for s in SOURCES:
        kinds = sorted({k for k, _, _ in s["files"]})
        n = len(s["files"])
        shown = s["files"][0][1] + (f" (+{n - 1} more)" if n > 1 else "")
        where = s["url"] or s["doi"]
        if s["doi"] and s["url"] and s["doi"] not in s["url"]:
            where = f"{s['url']}; DOI {s['doi']}"
        kind = " / ".join(KIND_LABEL[k] for k in kinds) + (" (rolling)" if s["rolling"] and "fetch" in kinds else "")
        out.append(f"| {s['short']} | {s['name']} | `{shown}` | {kind} | {where} | {s['accessed']} | "
                   f"{s['licence']} | {s['used_by']} |")
    out += ["", "| raw product (maintainer extractions only) | size | from | accessed | licence | note |",
            "|---|---|---|---|---|---|"]
    for r in RAW:
        mb = f"{r['mb']:g} MB" if r["mb"] else ""
        out.append(f"| `{r['rel']}` | {mb} | {r['url']} | {r['accessed']} | {r['licence']} | {r['note']} |")
    return "\n".join(out)




def render_compact() -> str:
    """The short table of the top-level README."""
    out = ["| source | kind | licence | feeds |", "|---|---|---|---|"]
    for s in SOURCES:
        kinds = sorted({k for k, _, _ in s["files"]})
        kind = " / ".join(KIND_LABEL[k] for k in kinds) + (" (rolling)" if s["rolling"] and "fetch" in kinds else "")
        out.append(f"| {s['short']}: {s['name']} | {kind} | {s['licence']} | {s['used_by']} |")
    return "\n".join(out)


# (file, begin marker, end marker, renderer)
README_BLOCKS = [
    (HERE / "README.md",
     "<!-- SOURCES TABLE: generated by `python tools/targets/target_inputs.py`; do not edit by hand -->",
     "<!-- END SOURCES TABLE -->", render_markdown),
    (REPO_ROOT / "README.md",
     "<!-- SOURCES TABLE (compact): generated by `python tools/targets/target_inputs.py --write-readme`; "
     "do not edit by hand -->", "<!-- END SOURCES TABLE (compact) -->", render_compact),
]


def readme_table_is_current() -> bool:
    for path, begin, end, render in README_BLOCKS:
        block = path.read_text().split(begin, 1)[1].split(end, 1)[0].strip("\n")
        if block != render():
            return False
    return True


def write_readme_table():
    for path, begin, end, render in README_BLOCKS:
        head, rest = path.read_text().split(begin, 1)
        _, tail = rest.split(end, 1)
        path.write_text(head + begin + "\n" + render() + "\n" + end + tail)


if __name__ == "__main__":
    import sys
    if "--write-readme" in sys.argv:
        write_readme_table()
    elif "--check-readme" in sys.argv:
        ok = readme_table_is_current()
        print("README sources table " + ("matches target_inputs.py" if ok else "is STALE: run with --write-readme"))
        sys.exit(0 if ok else 1)
    else:
        print(render_markdown())
