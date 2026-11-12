# data/sources/noaa — licence notes

| file | source | licence / conditions |
|---|---|---|
| `nasa_gmsl_annual.csv` | NOAA STAR Laboratory for Satellite Altimetry, global mean sea level, seasonal signals removed (`slr_sla_gbl_free_ref_90.csv`, release of 2026-05-20), reduced to calendar-year means (years with at least 20 records) by the research `download_obs.py`, as `tools/targets/fetch_inputs.py` does now | No licence. The LSA data-use policy (2008) requires this credit line: "Altimetry data are provided by the NOAA Laboratory for Satellite Altimetry". The products are experimental. |
| `noaa_thermosteric_w0-2000m_yearly.dat` | NOAA NCEI, World Ocean 0-2000 m thermosteric sea level, yearly (https://www.ncei.noaa.gov/access/global-ocean-heat-content/), accessed 2026-06-13 | No licence stated on the pages read (2026-11-12). Credit NOAA NCEI. |

Both are rolling products: NOAA replaces them with each release. These are the vintages the
calibration targets were built from; a download made today does not reproduce them (and so
does not reproduce `data/observations/calibration_targets.csv`). Shipped by Marcus's ruling of
2026-11-12. The file name `nasa_gmsl_annual.csv` is the research repository's; the data are
NOAA STAR's.
