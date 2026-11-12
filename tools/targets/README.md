# tools/targets: rebuild the calibration targets and drivers from their sources

Rebuilds the eight files in `data/` that come from observations, and proves the
rebuild matches the shipped files:

| shipped file | research file (SLR-RFF-BRICK) | built by |
|---|---|---|
| `data/observations/calibration_targets.csv` | `outputs/recalib_targets_ext.csv` | `prep_recalib_targets_ext.py` |
| `data/observations/calibration_targets_gsic_adj.csv` | `outputs/recalib_targets_ext_gsicadj.csv` | `ladrillo_data.py` |
| `data/drivers/glacier_block_temperature.csv` | `data/observations/t_glac_blocks.csv` | `ladrillo_data.py` |
| `data/constants/glacier_block_constants.csv` | `outputs/extc_block_constants.csv` | `ladrillo_data.py` |
| `data/drivers/greenland_zone_temperature.csv` | `data/observations/t_gis_zones.csv` | `build_t_gis.py` |
| `data/calibration/gis_amp_prior.csv` | `outputs/gis_amp_prior.csv` | `build_t_gis.py` |
| `data/constants/gis_amp_shape.csv` | `outputs/gis_amp_shape.csv` | `diag_gis_amp_cmip6.py` |
| `data/constants/gis_amp_shape_meta.csv` | `outputs/gis_amp_shape_meta.csv` | `diag_gis_amp_cmip6.py` |

The builders write the research file names; a last step applies the label-only
rename that made `data/` (`tools/import_l27_data.py`: SLOWP to SLOWG, FAST to
FASTG). Four intermediates the research repository also keeps are gated too:
`diag_constraint_anatomy_regions.csv`, `lws_grace_extension_L14.csv`,
`gis_driver_constants.csv`, `gis_amp_anchor.csv`.

**Status (2026-11-12, macOS arm64, Python 3.14.3, `requirements-targets.txt`).**
All eight shipped files and the four intermediates rebuild byte for byte (sha256
identical) from `data/sources/` plus the two Frederikse (2020) files fetched from
Zenodo, in about 6 s. Nothing else has to be fetched. One cell differs by design:
the `commit` column of `gis_amp_shape_meta.csv` records the commit of the code
that built the shape (`d59ce12` in the shipped file). With `--from-raw` the chain
starts from the raw gridded products instead, re-makes all seven extractions
byte for byte, and then rebuilds the same files (about 1-3 min). Gate: 20/20 and
27/27.

## Quick start

```bash
python -m venv .venv-targets && . .venv-targets/bin/activate
pip install -r tools/targets/requirements-targets.txt

python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs             # the plan: what to fetch, from where
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --download  # the two Frederikse (2020) files, Zenodo
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --verify

python tools/targets/build_targets.py --inputs-dir ~/ladrillo_inputs --work-dir /tmp/targets_build
python tools/targets/check_targets_identity.py --work-dir /tmp/targets_build
```

In the maintainer's research checkout every input is already present:

```bash
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --from-research <SLR-RFF-BRICK checkout>
python tools/targets/fetch_inputs.py --inputs-dir ~/ladrillo_inputs --extract
python tools/targets/build_targets.py --inputs-dir ~/ladrillo_inputs --work-dir /tmp/tb            # from the extractions
python tools/targets/build_targets.py --inputs-dir ~/ladrillo_inputs --work-dir /tmp/tb2 --from-raw # from the raw products
python tools/targets/check_targets_identity.py --work-dir /tmp/tb2
```

`build_targets.py` checks every input against `inputs.sha256` before any step
runs. One missing or mismatched file aborts the build, and the message names it.
It runs the scripts from a frozen copy (`<work-dir>/tools_snapshot/`) with BLAS
pinned to one thread, writes only under `--work-dir`, and records in
`<work-dir>/BUILD_PROVENANCE.json` every command and its wall time, the sha256 of
every script, input and output, the package versions, the Ladrillo commit, and
the one seed the chain uses (below). It touches `data/` only with `--install`.
Install only on purpose: the shipped posterior was calibrated on these targets.

### Rolling products and current vintages

Several sources are overwritten by their publishers (marked rolling below):
NOAA STAR, NOAA NCEI, Berkeley Earth, GISTEMP, HadCRUT5, the GRACE mascons and
CMIP6 on Pangeo. The vintages the targets were built from are no longer served,
so Ladrillo ships the small files made from them (ruling, 2026-11-12): the NOAA
STAR annual means and the NCEI thermosteric series as received, and extractions
of the others. The standard rebuild uses these and needs no archive.

To build on today's versions instead, `fetch_inputs.py` describes the route
(CURRENT VINTAGES in its plan): `--download` also writes the current STAR file,
reduced to annual means as the research did; place a current NCEI file the same
way; put current Berkeley Earth and GISTEMP grids in `raw/` and run `--extract
--no-pin-check`; then `build_targets.py --current-vintages`. Those files replace
the shipped copies unpinned (the build lists them). Today's vintages do not
reproduce the shipped bytes, and the gate reports which outputs move and by how
much.

## The chain

| # | step | script | research source | reads | writes (under `--work-dir`) |
|---|---|---|---|---|---|
| 1 | regions | `build_region_inventory.py` | `diag_constraint_anatomy.py` part 2 @ 3f7f2aa | GlacierMIP3 tables | `outputs/diag_constraint_anatomy_regions.csv` |
| 2 | lws | `build_lws_grace_extension.py` | `build_lws_grace_extension.py` @ d2e23d3 (time-series half) | GRACE land-water extraction, GlaMBIE, Frederikse workbook | `outputs/lws_grace_extension_L14.csv` |
| 3 | targets | `prep_recalib_targets_ext.py` | same @ bc50130, no flags | Frederikse, Dangendorf, GRACE, GlaMBIE, NOAA, step 2 | `outputs/recalib_targets_ext.csv` |
| 4 | glacier | `ladrillo_data.py` | same @ 7e142ff | step 3 (gsic), step 1, HadCRUT5 glacier extraction, GlaMBIE, GlacierMIP3 | `t_glac_blocks.csv`, `extc_block_constants.csv`, `recalib_targets_ext_gsicadj.csv` |
| 5 | greenland | `build_t_gis.py` | same @ d59ce12 (derivation half) | the three Greenland zone extractions, step 3 (gis) | `t_gis_zones.csv`, `gis_driver_constants.csv`, `gis_amp_prior.csv` |
| 6 | amp_anchor | `diag_gis_amp_anchor.py` | same @ 8b68d19 | step 5 | `outputs/gis_amp_anchor.csv` |
| 7 | amp_shape | `diag_gis_amp_cmip6.py` | same @ a8a1779, `--zone south` | 40 CMIP6 extractions, steps 5-6 | `gis_amp_shape.csv`, `gis_amp_shape_meta.csv` |

Maintainer extractions (`--from-raw`; `extract/`):

| extraction | script | research source | raw input |
|---|---|---|---|
| `hadcrut5/t_glac_hadcrut5.csv`, `t_glac_regions_hadcrut5.csv` | `extract/build_t_glac.py` | `build_t_glac.py` @ 597b7e1 | HadCRUT5, GTN-G regions, GlaMBIE areas |
| `greenland_zones/greenland_zone_series_<product>.csv` (3) | `extract/extract_greenland_zones.py` | `build_t_gis.py` @ d59ce12 (measurement half) | HadCRUT5, Berkeley Earth, GISTEMP, GTN-G region 05 |
| `grace/grace_lws_domain_monthly.csv` | `extract/extract_grace_lws_domain.py` | `build_lws_grace_extension.py` @ d2e23d3 (gridded half) | GRACE mascon grid, GTN-G regions |
| `dangendorf2024/dangendorf2024_gmsl_annual.csv` | `extract/extract_dangendorf_gmsl.py` | `diag_dangendorf_vs_frederikse.py` lines 83-94 @ a954701 | Dangendorf fields, global file, grid |
| `cmip6_gis/tas_series_gis_<model>.csv` (40) | `extract/reduce_cmip6_tas_gis.py` | same @ 8d52288 | CMIP6 on Pangeo (not re-run; see below) |
| `gmip3/gmip3_block_ladder_L27.csv` | not ported | `d1_multireservoir_cell.py` `block_ladder()` @ 903367f | GlacierMIP3 shifted-volume netCDF (1.47 GB) |

The two new extraction formats (the Greenland series and the GRACE monthly
series) are written at full round-trip precision, so the derivation steps
compute every number from the same doubles the research scripts used.

## Inputs and licences

`inputs.sha256` pins every input: `sources/` (shipped in `data/sources/`),
`fetch/` (fetch yourself into `--inputs-dir`) and `raw/` (raw products, read only
by the maintainer extractions). The table below is generated from
`target_inputs.py` (`python tools/targets/target_inputs.py`). Licences were read
from the publishers' pages and the Zenodo API on 2026-11-12; nothing was
downloaded.

Rulings (Marcus, 2026-11-12). A dataset whose licence leaves redistribution in
doubt is not shipped; it is fetched, and pinned by sha256. That applies to
Frederikse (2020): its original Zenodo record (3862995) is "Other (Not Open)",
and Dangendorf's CC BY 4.0 record (10621070), which redistributes the two files,
cannot relicense them. The four small derived files of rolling products are
shipped with licence notes, because their vintages are no longer served: NOAA
STAR (no licence; the data-use policy requires the credit line "Altimetry data
are provided by the NOAA Laboratory for Satellite Altimetry"), NOAA NCEI (no
licence stated), Berkeley Earth (CC BY-NC 4.0: non-commercial use only, not
covered by Ladrillo's MIT licence) and GISTEMP (no licence stated). Each has a
README next to it (`data/sources/noaa/`, `data/sources/greenland_zones/`). The
Dangendorf v2 file is shipped by ruling.

<!-- SOURCES TABLE: generated by `python tools/targets/target_inputs.py`; do not edit by hand -->
| source | what | files | kind | where from | accessed | licence (read 2026-11-12) | feeds |
|---|---|---|---|---|---|---|---|
| Frederikse et al. (2020) | sea-level budget: global basin time series (xlsx) and the 5000-member weighted ensemble (nc) | `frederikse2020/frederikse2020_global_basin_timeseries.xlsx (+1 more)` | fetch yourself | https://zenodo.org/records/10621070/files/global_basin_timeseries.xlsx?download=1 and https://zenodo.org/records/10621070/files/GMSL_ensembles_F20.nc?download=1 (rename on save); DOI 10.5281/zenodo.10621070 (redistributed by Dangendorf et al.; original 10.5281/zenodo.3862995) | 2026-05-21 (xlsx), 2026-07-19 (nc) | original record 3862995: 'Other (Not Open)', all rights reserved; record 10621070: CC BY 4.0, which cannot relicense Frederikse's files | every component target (ais, gis, gsic, steric, lws) 1900-2018 and its sigma; the total's budget-closure sigma; the land-water splice offset |
| Dangendorf et al. (2024), corrected global file | Kalman-smoother GMSL reconstruction, global series with standard errors (v2) | `dangendorf2024/dangendorf2024_KalmanSmootherHR_Global_v2.nc` | shipped | none: personal communication, S. Dangendorf to M. Sarofim, 2026-08-07 (fixes the slot shift in Zenodo 10621070's KalmanSmootherHR_Global.nc) | 2026-08-07 | none stated; shipped by Marcus's ruling (2026-11-12) | dang_sig, dang_lo, dang_hi 1900-2021 (the total target's own SE) |
| Dangendorf et al. (2024) | annual GMSL = cos(lat)-weighted mean of the HR field (extraction) | `dangendorf2024/dangendorf2024_gmsl_annual.csv` | shipped | https://zenodo.org/records/10621070; DOI 10.5281/zenodo.10621070 | 2026-07-20 | CC BY 4.0 (Zenodo API) | dang (the total target) 1900-2021 |
| JPL GRACE/GRACE-FO mascons RL06.3Mv04 | Antarctic and Greenland mass time series | `grace/grace_antarctica_mass.txt (+1 more)` | shipped (rolling) | https://podaac.jpl.nasa.gov/ (Earthdata login); DOI 10.5067/TEMSC-AT634, 10.5067/TEMSC-GT634 | 2026-06-13 | NASA EOSDIS data-use guidance: openly shared without restriction (CC0 unless marked) | ais and gis 2019-2026; the mask gates of the land-water extraction |
| JPL GRACE/GRACE-FO mascons RL06.3Mv04 | monthly mass over global land minus both ice sheets (extraction) | `grace/grace_lws_domain_monthly.csv` | shipped (rolling) | https://podaac.jpl.nasa.gov/ (Earthdata login); DOI 10.5067/TEMSC-3JC634 | 2026-08-24 | as above (NASA, without restriction) | lws 2019-2023, held flat to 2026 |
| GlaMBIE Team (2025) | regional calendar-year glacier mass change (results, calendar_years/ of glambie_data.zip) | `glambie/calendar_years/0_global.csv (+19 more)` | shipped | https://doi.org/10.5904/wgms-glambie-2024-07 | 2026-06-13 | WGMS data policy: open access with citation (= CC BY 4.0) | gsic 2019-2023 (global minus region 5); glacier block weights, modern rates, 2020 losses; the region-19 seam; glaciers removed from the land-water extension |
| NOAA NCEI | World Ocean 0-2000 m thermosteric sea level, yearly | `noaa/noaa_thermosteric_w0-2000m_yearly.dat` | shipped (rolling) | https://www.ncei.noaa.gov/access/global-ocean-heat-content/ (yearly thermosteric sea level, World Ocean, 0-2000 m) | 2026-06-13 | none stated on the pages read; shipped by ruling (rolling product, 2026-06 vintage) | steric 2019-2025 |
| NOAA STAR Laboratory for Satellite Altimetry | global mean sea level, seasonal signals removed; annual means of the 2026-05 release | `noaa/nasa_gmsl_annual.csv` | shipped (rolling) | https://www.star.nesdis.noaa.gov/socd/lsa/SeaLevelRise/slr/slr_sla_gbl_free_ref_90.csv (fetch_inputs.py --download derives the annual means as the research did) | 2026-05-20 | no licence; LSA data-use policy (2008) requires the credit line "Altimetry data are provided by the NOAA Laboratory for Satellite Altimetry"; shipped by ruling (rolling, 2026-05-20 release) | dang 2022-2024 |
| GlacierMIP3 (Zekollari et al. 2025) | Tables S1a and S3, region characteristics | `gmip3/table_S1a.csv (+2 more)` | shipped | https://zenodo.org/records/15046588; DOI 10.5281/zenodo.15046588 | 2026-08-07 | CC BY 4.0 (Zenodo API) | glacier inventory, response times, ISIMIP3 amplification |
| GlacierMIP3 (Zekollari et al. 2025) | committed-loss ladder per Ladrillo glacier reservoir (extraction) | `gmip3/gmip3_block_ladder_L27.csv` | shipped | https://zenodo.org/records/15046588; DOI 10.5281/zenodo.15046588 | 2026-08-07 | CC BY 4.0 (Zenodo API) | the glacier equilibrium curves (b, T_off) and their committed rungs |
| Met Office HadCRUT5.0.2.0 (Morice et al. 2021) | glacier-area-weighted annual temperature per RGI region, and the global mean (extraction) | `hadcrut5/t_glac_hadcrut5.csv (+1 more)` | shipped (rolling) | https://www.metoffice.gov.uk/hadobs/hadcrut5/; DOI 10.1029/2019JD032361 | 2026-08-06 | Open Government Licence v3 (Crown copyright); acknowledgement line required | the three glacier reservoir drivers; the observational amplification fit |
| Met Office HadCRUT5.0.2.0 | Greenland land-zone and global annual temperature (extraction) | `greenland_zones/greenland_zone_series_HadCRUT5.csv` | shipped (rolling) | https://www.metoffice.gov.uk/hadobs/hadcrut5/; DOI 10.1029/2019JD032361 | 2026-08-06 | Open Government Licence v3 (the zone weights use the Berkeley Earth 1-degree land fraction and the GTN-G region-05 polygon; no Berkeley values are in the file) | the Greenland driver; the amplification prior and anchor |
| Berkeley Earth | Greenland land-zone and global annual temperature (extraction) | `greenland_zones/greenland_zone_series_BerkeleyEarth.csv` | shipped (rolling) | https://berkeleyearth.org/data/ (Land_and_Ocean_LatLong1.nc) | 2026-08-09 | CC BY-NC 4.0 (berkeleyearth.org/data): non-commercial use only, not covered by Ladrillo's MIT licence; shipped by ruling (rolling, 2026-08 vintage) | the amplification prior and anchor (one of three products) |
| NASA GISS GISTEMP v4 | Greenland land-zone and global annual temperature (extraction) | `greenland_zones/greenland_zone_series_GISTEMP.csv` | shipped (rolling) | https://data.giss.nasa.gov/gistemp/ (gistemp1200_GHCNv4_ERSSTv5.nc.gz) | 2026-08-09 | none stated (citation of Lenssen et al. 2024 requested); shipped by ruling (rolling, 2026-08 vintage) | the amplification prior and anchor (one of three products) |
| CMIP6 (40 models) | Amon tas, annual global and Greenland-zone means, historical + ssp126/245/585 (extraction) | `cmip6_gis/tas_series_gis_ACCESS-CM2.csv (+39 more)` | shipped (rolling) | https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv (Pangeo) | 2026-08-12 / 2026-08-18 | CC BY 4.0 for all groups since 2022-10 (CMIP6 Terms of Use v6.2); GISS sources CC0 | the Greenland amplification shape S(dT) |

| raw product (maintainer extractions only) | size | from | accessed | licence | note |
|---|---|---|---|---|---|
| `HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc` | 32.5 MB | https://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.0.2.0/analysis/HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc | 2026-08-06 | OGL v3 | updated monthly in place |
| `Land_and_Ocean_LatLong1.nc` | 454.6 MB | https://berkeleyearth.org/data/ (Gridded, Land + Ocean, 1x1 degree) | 2026-08-09 | CC BY-NC 4.0 | updated in place |
| `gistemp1200_GHCNv4_ERSSTv5.nc` | 57 MB | https://data.giss.nasa.gov/pub/gistemp/gistemp1200_GHCNv4_ERSSTv5.nc.gz (gunzip) | 2026-08-09 | none stated | updated monthly in place |
| `GRCTellus.JPL.200204_202606.GLO.RL06.3M.MSCNv04CRI.nc` | 45.8 MB | https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4 (Earthdata login) | 2026-08-24 | NASA, without restriction | the name carries the last month; replaced monthly |
| `GlacReg_2023/GTN-G_202307_o1regions.shp` |  | https://doi.org/10.5904/gtng-glacreg-2023-07 (GlacReg_2023.zip) | 2026-08-06 | GTN-G data policy: open access with citation (= CC BY 4.0), also 'public domain for non-commercial use' |  |
| `GlacReg_2023/GTN-G_202307_o1regions.shx` |  | (same zip) | 2026-08-06 | as above |  |
| `GlacReg_2023/GTN-G_202307_o1regions.dbf` |  | (same zip) | 2026-08-06 | as above |  |
| `dangendorf2024_KalmanSmootherHR_Fields.nc` | 729.5 MB | https://zenodo.org/records/10621070/files/KalmanSmootherHR_Fields.nc?download=1 | 2026-07-20 | CC BY 4.0 | matches Zenodo's md5 |
| `dangendorf2024_KalmanSmootherHR_Global.nc` | 0.06 MB | https://zenodo.org/records/10621070/files/KalmanSmootherHR_Global.nc?download=1 | 2026-07-19 | CC BY 4.0 | matches Zenodo's md5; slot-shifted upstream, read only for the indicative-SE column |
| `dangendorf2024_grid_lonlat.csv` | 1.2 MB | (lon/lat of the Fields file's grid points; how it was produced is not recorded in the research repo) | 2026-07-20 | CC BY 4.0 (derived from Zenodo 10621070) | origin not recorded |
<!-- END SOURCES TABLE -->

The Frederikse files fetched from Zenodo 10621070 match the copies the targets
were built from (Zenodo's published md5, checked through its API on
2026-11-12). No other fetch was tested: `fetch_inputs.py --download` and
`--verify` are the first real test for each.

The HadCRUT5 and GISTEMP Greenland series and the 40 CMIP6 series are zone means
whose cell weights use the Berkeley Earth 1-degree land fraction and the GTN-G
region-05 polygon. No Berkeley Earth temperature is in them.

### Attribution

- HadCRUT5 (OGL v3), the acknowledgement the Met Office requires: "HadCRUT.5.0.2.0
  data were obtained from http://www.metoffice.gov.uk/hadobs/hadcrut5 on 6 August
  2026 and are © British Crown Copyright, Met Office 2026, provided under an Open
  Government License, http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/".
  Morice, C. P., et al. (2021), J. Geophys. Res. Atmos., 126, e2019JD032361.
- GlaMBIE Team (2025), Nature; data 10.5904/wgms-glambie-2024-07 (WGMS, CC BY 4.0).
- Zekollari, H., Schuster, L., et al. (2025), Science, 10.1126/science.adu4675;
  data 10.5281/zenodo.15046588 (CC BY 4.0).
- Dangendorf, S., et al. (2024); data 10.5281/zenodo.10621070 (CC BY 4.0). The
  v2 global file: S. Dangendorf, personal communication, 2026-08-07.
- NASA JPL GRACE/GRACE-FO mascons RL06.3Mv04 CRI, 10.5067/TEMSC-3JC634 (and the
  Antarctic and Greenland mass series, 10.5067/TEMSC-AT634, 10.5067/TEMSC-GT634).
- CMIP6: the modelling groups of the 40 models in `data/sources/cmip6_gis/`, under
  the CMIP6 Terms of Use (CC BY 4.0; GISS models CC0).
- NOAA STAR, the credit line its data-use policy requires: "Altimetry data are
  provided by the NOAA Laboratory for Satellite Altimetry".
- Berkeley Earth (CC BY-NC 4.0, non-commercial use only): Rohde, R. A., and
  Hausfather, Z. (2020), Earth Syst. Sci. Data 12, 3469-3479.
- GISTEMP v4: GISTEMP Team, GISS Surface Temperature Analysis, NASA GISS, and
  Lenssen, N., et al. (2024).
- NOAA NCEI World Ocean thermosteric sea level (no licence stated).
- Fetched, not shipped: Frederikse, T., et al. (2020), Nature 584, 393-397.

## What the port changed

Each script is a port, not a rewrite: the numerical code is copied verbatim and
each docstring has a PORT paragraph listing its changes. Across the chain:

- Hard-coded research paths (`prep_recalib_targets_ext.py:88`, `build_t_glac.py:45`)
  became `--sources` / `--work` / `--raw-dir` flags.
- `diag_constraint_anatomy.py` built the region table after `exec`-ing part of
  `d0_glacier_shootout.py`. Only its part 2 is ported; it reads three GlacierMIP3
  tables and nothing from the exec.
- The silent land-water fallback (`prep_recalib_targets_ext.py:483-497`: a missing
  GRACE extension meant LWS held flat from 2018, a different target, with a
  warning) is now an error.
- The circular read is gone. The land-water builder took the Frederikse TWS for
  its splice offset from `recalib_targets_ext.csv`, which the targets builder
  writes after reading the land-water output. It now computes the TWS from the
  Frederikse workbook with the targets builder's own lines. One detail of the old
  path is kept: the values go through an in-memory CSV write and read with
  pandas' defaults, because pandas' default parser shifts some 17-digit decimals
  by one ulp (5.6e-17 here). Without that, the extension moves by up to 1.1e-16
  cm and the targets are no longer byte-identical (mutant M2 below).
- The GRACE mascon file is named and pinned. The research script globbed for the
  newest `GRCTellus.JPL.*` file, so its output changed whenever a newer month was
  downloaded next to it.
- `build_t_gis.py` compared its rebuilt driver with its own previous output
  before overwriting it (lines 362-377): a gate reading the path its script
  writes. The identity gate replaces it.
- The research `diag_gis_amp_cmip6.py` raised a NameError while drawing its figure
  (`lf` is never imported), after the shape files were written. The figure is not
  ported.
- GlaMBIE files are read from `data/sources/glambie/calendar_years/` instead of out
  of `glambie_data.zip`: the same members, byte for byte.

### Dropped or archived

Not ported because nothing Ladrillo reads depends on them: the IMBIE-2026 AIS
target branch (`--ais-imbie2026`, the closed L28-L35 arms) and the IMBIE 2023
cross-check with its two input files; the splice diagnostic PNG and the
`recalib_targets_ext_sources.csv` sidecar; `t_gis_zones_allproducts.csv`; the
provenance markdown files, console tables and figures of `build_t_glac.py`,
`build_t_gis.py` and `build_lws_grace_extension.py`; the per-window table, the
bootstrap curves, the summary, the `gis_amp_shape_fullcurve` arm and the other
zones of `diag_gis_amp_cmip6.py` (their bootstraps, seed 2026, never reached the
shape, a PCHIP through bin medians); part 1 of `diag_constraint_anatomy.py`. The
research scripts stay in SLR-RFF-BRICK as the record.

### Not ported: built once, frozen

- `gmip3/gmip3_block_ladder_L27.csv`: the 12 rows of the research ladder cache
  (`outputs/d1_block_ladder_cache.csv` @ 3165dca) for the three region sets the
  glacier reservoirs use. The cache was filled by `d1_multireservoir_cell.py` from
  the GlacierMIP3 shifted-volume netCDF (1.47 GB) with quantile LOWESS (`moepy`);
  that code is part of a development cell, not ported. The other rows (region sets
  tried and dropped in 2026-08) are not shipped.
- `cmip6_gis/` (40 files): `extract/reduce_cmip6_tas_gis.py` is kept as the
  record but not re-run. It streams from the public Pangeo archive, whose stores
  are added, retracted and re-versioned, so a re-run is not expected to give the
  same bytes.
- No builder exists in any history for four constants, shipped as frozen data
  (`data/MANIFEST.csv` gives their research source and commit):
  `constants/central_row.csv` (the medoid draw, `medoid_post_idx`, of the June 2026
  BRICK-Mengel posterior; research `outputs/recalib_central_row.csv` @ f2b0a8d);
  `constants/paleo_dais_marginals.csv` and `constants/paleo_fastdyn_draws.csv`
  (from the MimiBRICK DAIS fast-dynamics paleo ensemble; @ 20aaabf and 15d84e5);
  `constants/paleo_geo_prior_ton.csv` (@ a954701; the script that made it,
  `compute_paleo_geo_prior_ton.jl`, is in no git history). The calibration's start
  rows and proposal covariance (`data/calibration/`) come from the L26 run and
  cannot be rebuilt without it.

## Seeds

One step draws random numbers: `ladrillo_data.py` seeds each glacier
reservoir's multi-start Nelder-Mead jitter with `FIT_RNG_SEED = 2026`, one
stream per (block, inventory basis, amplification) via `block_rng()`. It affects
the fitted `b`, `T_off`, `kappa` and `nu` columns of `glacier_block_constants.csv`
(mutant M3 moves them). The output files are byte-gated research files and
cannot carry a new column, so the seed is recorded in `BUILD_PROVENANCE.json`.

## The identity gate

`check_targets_identity.py --work-dir DIR` compares by sha256, with no tolerance:

- **shipped**: each of the eight `data/` files against `data/MANIFEST.csv`'s
  sha256, and its research-named twin against MANIFEST's `source_sha256`. Each
  reference `data/` file is first checked against MANIFEST, so a reference edited
  in place cannot pass.
- **intermediate**: the four research intermediates against `reference.sha256`.
- **extraction** (after `--from-raw`): the seven re-made extractions against
  `inputs.sha256` and the shipped copies.

The one exemption is the `commit` cell of `gis_amp_shape_meta.csv` (verdict
PASS-LABEL): every other byte of that file must be identical. On a failure the
gate reports whether the values are equal (a formatting-only difference still
fails) or the largest |difference| and where. `--self-test` checks the
comparison: a one-ulp change, a dropped row, CRLF line endings, a changed meta
value, and a changed meta value next to a changed label all FAIL; a changed label
alone is PASS-LABEL.

`mutation_test_targets.py --inputs-dir DIR` rebuilds with one change at a time
and runs the gate. Result of 2026-11-12 (`mutation_result.txt`):

```
targets gate mutation test, 2026-11-12 09:13, Ladrillo ab3d1e8 + the 2026-11-12d working tree (run from a frozen copy), maintainer inputs (fetch_inputs.py --from-research)
ok  M1 altimetry sigma 4.0 -> 4.1 mm                               expect FAIL      got FAIL        TARGETS IDENTITY: 18/20 pass, 2 FAIL (2 on the named provenance label); failing: data/observations/calibration_targets.csv, outputs/recalib_targets_ext.csv
ok  M2 land-water read-back emulation removed (1-ulp)              expect FAIL      got FAIL        TARGETS IDENTITY: 17/20 pass, 3 FAIL (2 on the named provenance label); failing: data/observations/calibration_targets.csv, outputs/recalib_targets_ext.csv, outputs/lws_grace_extension_L14.csv
ok  M3 glacier fit seed 2026 -> 2027                               expect FAIL      got FAIL        TARGETS IDENTITY: 18/20 pass, 2 FAIL (2 on the named provenance label); failing: data/constants/glacier_block_constants.csv, outputs/extc_block_constants.csv
ok  M4 amp shape support to 3.25 K (the fullcurve arm)             expect FAIL      got FAIL        TARGETS IDENTITY: 16/20 pass, 4 FAIL (0 on the named provenance label); failing: data/constants/gis_amp_shape.csv, outputs/gis_amp_shape.csv, data/constants/gis_amp_shape_meta.csv, outputs/gis_amp_shape_meta.csv
ok  M5 Greenland driver from Berkeley Earth instead of HadCRUT5    expect FAIL      got FAIL        TARGETS IDENTITY: 18/20 pass, 2 FAIL (2 on the named provenance label); failing: data/drivers/greenland_zone_temperature.csv, data/observations/t_gis_zones.csv
ok  C1 control: a console message only                             expect PASS      got PASS        TARGETS IDENTITY: 20/20 pass (2 on the named provenance label); failing: -
ok  C2 control: melt-rate smoothing 11 -> 9 yr (feeds only an intermediate) expect FAIL-ONLY got FAIL        TARGETS IDENTITY: 19/20 pass, 1 FAIL (2 on the named provenance label); failing: outputs/gis_driver_constants.csv
ok  I1 one digit in a shipped input (GlaMBIE global)               expect ABORT     got ABORT       ABORT: inputs do not match inputs.sha256 (fetch_inputs.py --verify says more); nothing was run
ok  I2 one byte in a fetched input (Frederikse workbook)           expect ABORT     got ABORT       ABORT: inputs do not match inputs.sha256 (fetch_inputs.py --verify says more); nothing was run
ok  I3 one digit in a shipped rolling-product file (NOAA STAR annual means) expect ABORT     got ABORT       ABORT: inputs do not match inputs.sha256 (fetch_inputs.py --verify says more); nothing was run
MUTATION TEST PASS: every mutant behaved as expected
```

## Runtime

Apple M4, 10 cores: about 6 s from the extractions (no step takes over 2 s);
about 50 s with `--from-raw` at a load average near 3, and 2-3 min under load
(the Greenland and GRACE extractions dominate, and sha256 of the 1.4 GB of raw
inputs). Memory stays under 2 GB.
