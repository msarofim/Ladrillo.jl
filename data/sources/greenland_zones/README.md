# data/sources/greenland_zones — licence notes

Annual Greenland land-zone (south, all, central, north) and global-mean temperature anomalies,
one file per product, made by `tools/targets/extract/extract_greenland_zones.py` from the gridded
products. Zone weights: GTN-G 2023 region 05 polygon x the Berkeley Earth 1-degree land fraction.

| file | product (vintage) | licence / conditions |
|---|---|---|
| `greenland_zone_series_HadCRUT5.csv` | HadCRUT.5.0.2.0 analysis ensemble mean (downloaded 2026-08-06) | Open Government Licence v3, © British Crown Copyright, Met Office; the acknowledgement in `tools/targets/README.md` |
| `greenland_zone_series_BerkeleyEarth.csv` | Berkeley Earth Land and Ocean 1-degree grid (downloaded 2026-08-09) | **CC BY-NC 4.0** (berkeleyearth.org/data): attribution, non-commercial use only. This file is NOT covered by Ladrillo's MIT licence. Cite Rohde and Hausfather (2020), Earth Syst. Sci. Data 12, 3469-3479. |
| `greenland_zone_series_GISTEMP.csv` | NASA GISS GISTEMP v4, 1200 km (downloaded 2026-08-09) | No licence stated (2026-11-12); GISS asks for citation of the GISTEMP page and Lenssen et al. (2024). |

All three products are updated in place by their publishers; a download made today does not
reproduce these files. The Berkeley Earth and GISTEMP files are shipped by Marcus's ruling of
2026-11-12. All three feed the Greenland amplification prior (`data/calibration/gis_amp_prior.csv`)
and its anchor; the HadCRUT5 file also feeds the Greenland driver.
