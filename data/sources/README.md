# data/sources

The source files of Ladrillo's calibration targets and drivers that Ladrillo.jl ships:
the inputs of `tools/targets/` (see `tools/targets/README.md` for the chain, the gate,
and the full sources table with licences).

Shipped here only when the licence clearly allows redistribution (attribution at
most). Files whose licence leaves redistribution in doubt are not here; fetch them
with `tools/targets/fetch_inputs.py` (Frederikse 2020, NOAA STAR, NOAA NCEI, and the
Berkeley Earth and GISTEMP Greenland series).

| directory | what | licence |
|---|---|---|
| `dangendorf2024/` | the corrected global file (v2, personal communication 2026-08-07) and the annual GMSL extraction of the Zenodo 10621070 fields | v2: none stated, shipped by ruling; extraction: CC BY 4.0 |
| `grace/` | JPL GRACE/GRACE-FO mascon RL06.3Mv04: Antarctic and Greenland mass series, and the monthly land-water extraction | NASA, without restriction |
| `glambie/calendar_years/` | GlaMBIE 2025 regional calendar-year glacier mass change (members of glambie_data.zip) | CC BY 4.0 (WGMS) |
| `gmip3/` | GlacierMIP3 Tables S1a, S3, region characteristics, and the committed-loss ladder for Ladrillo's three glacier reservoirs | CC BY 4.0 |
| `hadcrut5/` | HadCRUT5.0.2.0 glacier-area-weighted temperature per RGI region and global (extraction) | OGL v3, Crown copyright |
| `greenland_zones/` | HadCRUT5.0.2.0 Greenland land-zone and global temperature (extraction) | OGL v3, Crown copyright |
| `cmip6_gis/` | 40 CMIP6 models: annual global and Greenland-zone tas (extraction) | CC BY 4.0 (GISS: CC0) |

`tools/targets/inputs.sha256` pins every file. Attribution lines are in
`tools/targets/README.md`.
