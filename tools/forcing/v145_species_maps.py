"""
v145_species_maps.py -- the three species maps build_emissions_v160_ssp245.py
needs, VENDORED VERBATIM from the research repo's
FaIRtoFrEDI/build_emissions_v145_ssp245.py @ 8984a13.

The research builder parsed these out of that script with `ast` so the two
builders could not drift apart. That script is private and frozen provenance,
so the definitions are copied here as source text (extracted with
ast.get_source_segment, not retyped). `build_emissions_v160_ssp245.py
--v145-builder PATH` re-parses the original and refuses to run if this copy has
drifted from it.
"""

AR6_TO_V145 = {
    "FossilCO2":  "CO2 FFI",
    "OtherCO2":   "CO2 AFOLU",
    "CH4":        "CH4",
    "N2O":        "N2O",
    "SOx":        "Sulfur",
    "CO":         "CO",
    "NMVOC":      "VOC",
    "NOx":        "NOx",
    "BC":         "BC",
    "OC":         "OC",
    "NH3":        "NH3",
    "CF4":        "CF4",
    "C2F6":       "C2F6",
    "C6F14":      "C6F14",
    "HFC23":      "HFC-23",
    "HFC32":      "HFC-32",
    "HFC43_10":   "HFC-4310mee",
    "HFC125":     "HFC-125",
    "HFC134a":    "HFC-134a",
    "HFC143a":    "HFC-143a",
    "HFC227ea":   "HFC-227ea",
    "HFC245fa":   "HFC-245fa",
    "SF6":        "SF6",
    "CFC_11":     "CFC-11",
    "CFC_12":     "CFC-12",
    "CFC_113":    "CFC-113",
    "CFC_114":    "CFC-114",
    "CFC_115":    "CFC-115",
    "CARB_TET":   "CCl4",
    "MCF":        "CH3CCl3",
    "HCFC_22":    "HCFC-22",
    "HCFC_141B":  "HCFC-141b",
    "HCFC_142B":  "HCFC-142b",
    "HALON1211":  "Halon-1211",
    "HALON1202":  None,    # not in v1.4.5; drop
    "HALON1301":  "Halon-1301",
    "HALON2402":  "Halon-2402",
    "CH3BR":      "CH3Br",
    "CH3CL":      "CH3Cl",
}

GROUP_A_TEMPLATE_AR6 = "FossilCO2"

V145_ORPHAN_SPECIES = {
    "C3F8", "C4F10", "C5F12", "C7F16", "C8F18", "c-C4F8",
    "HFC-152a", "HFC-236fa", "HFC-365mfc",
    "NF3", "SO2F2",
    "CH2Cl2", "CHCl3",
}
