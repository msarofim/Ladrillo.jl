# data/sources

The source files of Ladrillo's calibration targets and drivers that Ladrillo.jl ships:
the inputs of `tools/targets/` (see `tools/targets/README.md` for the chain, the gate,
and the full sources table with licences).

Every file keeps its source's licence; the table says which. The files are data, not
code: Ladrillo's MIT licence does not cover them, and the Berkeley Earth file is for
non-commercial use only (CC BY-NC 4.0). Directories whose files carry a condition have
their own README with the licence note (`noaa/`, `greenland_zones/`).

Not here: the two Frederikse (2020) files, whose original record reserves all rights.
Fetch them from Zenodo 10621070 with `tools/targets/fetch_inputs.py --download`.

The NOAA files and the Berkeley Earth and GISTEMP series come from rolling products
that their publishers overwrite. They are shipped (ruling, 2026-11-12) because the
vintages they were made from are no longer served; a download made today does not
reproduce them.

| directory | what | licence |
|---|---|---|
| `dangendorf2024/` | the corrected global file (v2, personal communication 2026-08-07) and the annual GMSL extraction of the Zenodo 10621070 fields | v2: none stated, shipped by ruling; extraction: CC BY 4.0 |
| `grace/` | JPL GRACE/GRACE-FO mascon RL06.3Mv04: Antarctic and Greenland mass series, and the monthly land-water extraction | NASA, without restriction |
| `glambie/calendar_years/` | GlaMBIE 2025 regional calendar-year glacier mass change (members of glambie_data.zip) | CC BY 4.0 (WGMS) |
| `gmip3/` | GlacierMIP3 Tables S1a, S3, region characteristics, and the committed-loss ladder for Ladrillo's three glacier reservoirs | CC BY 4.0 |
| `hadcrut5/` | HadCRUT5.0.2.0 glacier-area-weighted temperature per RGI region and global (extraction) | OGL v3, Crown copyright |
| `greenland_zones/` | Greenland land-zone and global temperature, one file per product (extractions): HadCRUT5.0.2.0, Berkeley Earth, GISTEMP v4 | HadCRUT5: OGL v3, Crown copyright; Berkeley Earth: CC BY-NC 4.0 (non-commercial); GISTEMP: none stated (see `greenland_zones/README.md`) |
| `noaa/` | NOAA STAR altimetry GMSL, annual means of the 2026-05-20 release; NOAA NCEI 0-2000 m thermosteric sea level, yearly | STAR: no licence, credit line required ("Altimetry data are provided by the NOAA Laboratory for Satellite Altimetry"); NCEI: none stated (see `noaa/README.md`) |
| `cmip6_gis/` | 40 CMIP6 models: annual global and Greenland-zone tas (extraction) | CC BY 4.0 (GISS: CC0) |

`tools/targets/inputs.sha256` pins every file. Attribution lines are in
`tools/targets/README.md`.
