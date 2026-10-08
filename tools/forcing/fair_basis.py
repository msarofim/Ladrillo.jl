# PORT (Ladrillo tools/forcing/, 2026-10-08): copied VERBATIM from FaIRtoFrEDI
# scripts/fair_basis.py @ 8984a13. No code change below this line; on 2026-10-08
# comments and one printed path naming private notes / research paths were trimmed.
"""
fair_basis.py — declare, in ONE place, the two input choices a FaIR run makes
that cannot be read off its output: which EMISSIONS BASIS it is on, and which
part of its history is PRESCRIBED FORCING rather than raw emissions.

WHY. Two failure modes, both of which have already happened in this tree.

1. THE BASIS SWAP IS SILENT. The v1.6.0
   posterior is conditioned on CMIP7 historical emissions; pairing it with
   RCMIP-native or Smith-2024 history puts a -0.05 to -0.12 K cold bias in and
   raises no error. It cost a retracted result in Aug 2026 and, measured on the
   van Vuuren markers 2026-10-02, still moves the 2015-2024 GMST bias vs IGCC
   from -0.109 to -0.061 K. ⚠ It does NOT move an absolute-threshold statistic:
   P(ERF>8.5) @2150 on marker H is 42.212 % on both bases. A first pass here
   claimed 6.4 pp and was RETRACTED the same day -- it had re-referenced ERF to
   1850-1900 before comparing it to 8.5 W/m2, which double-references a series
   FaIR already reports relative to 1750. A filename infix records the basis; this module makes a
   driver REFUSE the incoherent pair instead of merely recording it.

2. "DID THAT RUN USE PRESCRIBED FORCING OR EMISSIONS?" IS NOT ANSWERABLE FROM
   THE OUTPUT. Every FaIR run here is MIXED: a handful of species arrive as a
   prescribed ERF time series and the rest are run through the gas cycle from
   emissions. Which handful depends on the calibration version AND on the
   calibration directory (calibration_v160 runs Land use `calculated`;
   calibration_v160_prod runs it prescribed per CMIP7 marker). Nothing in a
   `fair_mean_gmst_*.csv` says which. `declare()` returns the one-line
   declaration a driver prints and stamps into its provenance string.

NEITHER FUNCTION TOUCHES A FaIR OBJECT. They read the species properties and
the parameter table, so they can run before `allocate()` and fail at startup.

USAGE in a driver, after the calibration paths are resolved:

    from fair_basis import declare, require_basis

    basis = require_basis(emis_path, calib_version=CAL_VER)   # exits if incoherent
    decl  = declare(species, properties, forcing_file, df_configs)
    print(decl.text)
    provenance = f"{MODEL_LABEL}; {basis.text}; {decl.text}; ..."
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# -- Named constants ---------------------------------------------------------
# The provenance infix that asserts "CMIP7 historical verbatim + a future
# harmonized onto it at the anchor". Written by the build_emissions_v160_*
# builders; see any of their module docstrings for the construction.
CMIP7_BASIS_SLUG = "cmip7harm"
# `build_emissions_v160_ssp245.py` predates the slug and spells the same
# construction out longhand. Both are the CMIP7 basis.
CMIP7_BASIS_ALIASES = ("cmip7harm", "cmip7_ssp245_harmonized")
# Calibration versions whose posterior was conditioned on CMIP7 history, and
# which therefore MUST be driven with a CMIP7-basis emissions file.
CMIP7_CONDITIONED_CALIBS = ("1.6.0",)

HIST_ANCHOR_YEAR = 2023      # CMIP7 history ends 2023.5; the harmonization join

_LEGACY_BASES = {
    "smith2020": "Smith 2024 historical (calib 1.4.5 lineage)",
    "smith2024": "Smith 2024 historical (calib 1.4.5 lineage)",
    "spliced_ext": "Smith 2024 historical + van Vuuren future (calib 1.4.5 lineage)",
    "AR6_emissions": "AR6 / Rennels elemental-mass file (v1.4.1 lineage)",
    "_raw": "RCMIP-native, unharmonized",
}


@dataclass
class Basis:
    """What emissions a run is on, and whether that is coherent with the calib."""
    path: Path
    calib_version: str
    is_cmip7: bool
    label: str
    waived: bool = False

    @property
    def text(self) -> str:
        tag = "CMIP7-basis" if self.is_cmip7 else "LEGACY-BASIS"
        w = " [WAIVED]" if self.waived else ""
        return f"emissions {self.path.name} ({tag}: {self.label}){w}"


@dataclass
class ForcingDeclaration:
    """Which species are prescribed forcing and which are emissions-driven."""
    prescribed: list
    emissions_driven: list
    concentration_driven: list
    calculated: list
    forcing_file: Path | None
    scaled: list = field(default_factory=list)

    @property
    def text(self) -> str:
        pf = ",".join(self.prescribed) or "none"
        extra = ""
        if self.concentration_driven:
            extra = f"; concentration-driven {','.join(self.concentration_driven)}"
        src = self.forcing_file.name if self.forcing_file else "(no forcing file)"
        return (f"history basis MIXED: {len(self.emissions_driven)} species "
                f"emissions-driven through the gas cycle, "
                f"{len(self.prescribed)} prescribed as ERF from {src} "
                f"[{pf}], {len(self.calculated)} calculated{extra}")


def classify_basis(emis_path, calib_version) -> Basis:
    """Read the basis off the emissions FILENAME. No gate, no exit."""
    p = Path(emis_path)
    name = p.name
    if any(a in name for a in CMIP7_BASIS_ALIASES):
        return Basis(p, calib_version, True,
                     f"CMIP7 historical 1750.5-{HIST_ANCHOR_YEAR}.5 verbatim + "
                     f"a future harmonized at {HIST_ANCHOR_YEAR}.5")
    for key, lab in _LEGACY_BASES.items():
        if key in name or key in str(p.parent):
            return Basis(p, calib_version, False, lab)
    return Basis(p, calib_version, False, "UNRECOGNISED — the filename carries no "
                                          "basis marker, so the basis cannot be "
                                          "audited from the output")


def require_basis(emis_path, calib_version, allow_legacy=False) -> Basis:
    """Gate a (calibration, emissions) pair. Exits on the incoherent combination.

    A calibration in CMIP7_CONDITIONED_CALIBS driven with a non-CMIP7 emissions
    file is the documented input mismatch, not a variant — so the default is to
    refuse. `allow_legacy=True` (wire it to an explicit CLI flag, never to a
    default) permits it and marks the Basis WAIVED, so the waiver propagates
    into the provenance string rather than being lost at the call site.
    """
    basis = classify_basis(emis_path, calib_version)
    if str(calib_version) in CMIP7_CONDITIONED_CALIBS and not basis.is_cmip7:
        if not allow_legacy:
            sys.exit(
                f"ERROR: calib {calib_version} was conditioned on CMIP7 historical "
                f"emissions, but {basis.path.name} is {basis.label}.\n"
                f"  That pair is the measured input mismatch (-0.05 to -0.12 K of "
                f"historical GMST fit; it does NOT move an absolute-threshold "
                f"statistic), not a variant.\n"
                f"  Build a CMIP7-basis file -- tools/forcing/build_emissions_v160_*.py "
                f"write an '{CMIP7_BASIS_SLUG}' infix -- or pass the driver's "
                f"explicit legacy-basis flag and say so wherever the number is used."
            )
        basis.waived = True
        print(f"  ⚠ LEGACY EMISSIONS BASIS WAIVED: {basis.path.name} is "
              f"{basis.label} and calib {calib_version} expects the CMIP7 basis. "
              f"This run carries the input mismatch BY REQUEST.")
    return basis


def declare(species, properties, forcing_file=None, df_configs=None,
            forcing_array=None, species_index=None) -> ForcingDeclaration:
    """Declare the prescribed-vs-emissions split, and gate it if given the data.

    `properties` is `fair.io.read_properties()`'s second return. When
    `forcing_array` (timebounds x config x specie, as loaded) and
    `species_index` are supplied, a prescribed specie that reads back all-NaN is
    an error: production must not silently drop a forcing term. When
    `df_configs` is supplied, a prescribed specie with no `forcing_scale[...]`
    column is an error: FaIR applies forcing_scale internally to CALCULATED
    species only, so an unscaled prescribed term is indistinguishable from a
    cooler posterior.
    """
    def by_mode(m):
        return [s for s in species if properties[s]["input_mode"] == m]

    d = ForcingDeclaration(
        prescribed=by_mode("forcing"),
        emissions_driven=by_mode("emissions"),
        concentration_driven=by_mode("concentration"),
        calculated=by_mode("calculated"),
        forcing_file=Path(forcing_file) if forcing_file else None,
    )
    if d.prescribed and d.forcing_file is None:
        sys.exit(f"ERROR: {len(d.prescribed)} species are input_mode='forcing' "
                 f"({','.join(d.prescribed)}) but no forcing file was named.")
    if forcing_array is not None and species_index is not None:
        blank = [s for s in d.prescribed
                 if np.isnan(forcing_array[..., species_index[s]]).all()]
        if blank:
            sys.exit(f"ERROR: {d.forcing_file.name} has no rows for prescribed "
                     f"species {blank}; production cannot run with a missing "
                     f"forcing term.")
    if df_configs is not None:
        missing = [s for s in d.prescribed
                   if f"forcing_scale[{s}]" not in df_configs.columns]
        if missing:
            sys.exit(f"ERROR: prescribed species {missing} have no "
                     f"forcing_scale[...] column; they would run UNSCALED.")
        d.scaled = list(d.prescribed)
    return d
