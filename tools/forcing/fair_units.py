# PORT (Ladrillo tools/forcing/, 2026-10-08): copied VERBATIM from FaIRtoFrEDI
# scripts/fair_units.py @ 8984a13. No code change below this line; on 2026-10-08
# comments naming private research branches were trimmed.
"""
fair_units.py — write FaIR's OWN spelling of an emissions unit, never the source
file's, and refuse any relabel that is not provably factor-1.0.

WHY. CMIP7's historical emissions frame spells HFC-4310mee's unit
`kt HFC4310/yr`; FaIR, RCMIP and the v1.4.5 payloads all spell it
`kt HFC4310mee/yr`. Every v160 emissions builder carried the CMIP7 string
through to its output.

⚠ NO NUMBER IS WRONG, AND THE FAILURE MODE IS NOT WHAT IT LOOKS LIKE.
`fair.io.fill_from._parse_unit` treats an unrecognised compound as USER-DEFINED:
it registers `compound_convert["HFC4310"] = {"HFC4310": 1}`, so the factor is
exactly 1.0 and values load untouched. But it ALSO mutates the module-level
`desired_emissions_units["HFC-4310mee"]`, and that global stays mutated, so a
later read of a correctly-spelled source in the same interpreter dies:
    v160 then v145   -> KeyError: 'HFC4310'
    v160 then RCMIP  -> KeyError: 'HFC4310'
while v145-then-v160 is fine. An order-dependent KeyError in a shared global is
what surfaces months later in somebody else's script.

PROVENANCE OF THIS IMPLEMENTATION. The TWO-WARRANT check below is
`build_emissions_v160_ssp245.py`'s, written independently on one research
branch 2026-09-24 and promoted here verbatim in logic 2026-10-02. It is kept in
preference to the version on another research branch, which required only
that the mass prefix and the time unit match and then relabelled. That is
NECESSARY BUT NOT SUFFICIENT: `kt HFC32/yr` and `kt HFC125/yr` share prefix and
time, so the weaker gate would have silently relabelled one gas as another.

THE TWO WARRANTS, and why one alone cannot do the job:
  machine — FaIR's own `compound_convert` prices the pair at exactly 1.0. This
            covers NOx->NO2, but CANNOT cover the case this function exists for:
            FaIR has no `HFC4310` entry at all, that being the whole defect.
  human   — the pair is in UNIT_ALIASES, where the equivalence is ARGUED in a
            comment rather than looked up.
Anything in NEITHER is refused, not relabelled.

⚠ An alias set alone is NOT a fix. Both v160 SSP builders carried UNIT_ALIASES
and used it only to make a comparison GATE permissive, then wrote the
non-canonical spelling to disk anyway. A hand-maintained equivalence table
records a discrepancy; only normalising at the point of authorship resolves it.

The CMIP7 cache is NOT edited — it is a faithful copy of upstream, and
normalising belongs where we author a FaIR input.
"""

from __future__ import annotations

import sys

# Human-verified factor-1.0 label equivalences. Each is an ARGUED claim:
#   Mt NO2/yr  == Mt NOx/yr      CMIP7 reports NOx mass as NO2; FaIR's NOx
#                                species expects exactly that convention.
#   kt HFC4310/yr == kt HFC4310mee/yr   same gas (HFC-43-10mee); CMIP7 drops the
#                                `mee` suffix. FaIR cannot price this pair
#                                because it has no `HFC4310` entry.
UNIT_ALIASES = {
    ("Mt NO2/yr", "Mt NOx/yr"),
    ("kt HFC4310/yr", "kt HFC4310mee/yr"),
}


def fair_canonical_unit(specie, unit, aliases=UNIT_ALIASES):
    """FaIR's own spelling for `specie`, or `unit` where FaIR has no opinion."""
    from fair.structure.units import compound_convert, desired_emissions_units

    unit = str(unit).strip()
    want = desired_emissions_units.get(specie)
    if want is None or want == unit:
        return unit
    have_c = unit.split()[-1].split("/")[0]
    want_c = want.split()[-1].split("/")[0]
    machine_ok = compound_convert.get(have_c, {}).get(want_c) == 1.0
    human_ok = (unit, want) in aliases or (want, unit) in aliases
    if unit.split()[0] != want.split()[0] or not (machine_ok or human_ok):
        sys.exit(f"ERROR: {specie} is '{unit}' but FaIR wants '{want}', and that "
                 f"is NOT a verified factor-1.0 relabel (FaIR cannot price it and "
                 f"UNIT_ALIASES does not claim it). Convert the VALUES or keep the "
                 f"unit -- never relabel a quantity you have not rescaled.")
    return want


def report(renamed):
    """One line per specie whose unit string was normalised."""
    if not renamed:
        print("  [UNIT-CANON] no unit string needed normalising")
        return
    print(f"  [UNIT-CANON] normalised {len(renamed)} unit string(s) to FaIR's "
          f"spelling (values untouched):")
    for sp, (got, want) in sorted(renamed.items()):
        print(f"    {sp:<14} {got!r} -> {want!r}")
