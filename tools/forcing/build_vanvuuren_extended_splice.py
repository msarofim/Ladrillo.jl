"""
build_vanvuuren_extended_splice.py

Build van Vuuren spliced emissions CSVs that use the PUBLISHED ScenarioMIP
post-2100 EXTENSIONS (van Vuuren et al. 2026, GMD 19:2627-2656, Sect. 4 +
Table 3, Fig. 2) instead of the constant-at-2100 placeholder tail written by
`run_fair_vanvuuren_erf.py`.

WHY THIS EXISTS
---------------
`run_fair_vanvuuren_erf.py` holds every species constant at its 2100.5 value
through 2300.5. That was a placeholder: the IIASA `scenariomip-2026-prerelease`
pull covers 2023-2100 only. The paper DOES define extensions to 2500 for all
seven markers, and the final extended emissions are published on Zenodo. Any
post-2100 statement (e.g. P(ERF>8.5) at 2150) must use the real extension.

EXTENSION SOURCE
----------------
Zenodo 10.5281/zenodo.20713982 (Sanderson & Smith), record 20713982, version
**1.1.1**, published 2026-06-16, "Data used in the figures of the paper ...".
File `data/fair-inputs/emissions_1750-2500.csv` inside
`chrisroadmap/cmip7-scenariomip-1.1.1.zip`.

NOTE the version: the DOI printed in the paper's Code-and-data-availability
section (10.5281/zenodo.14382495) resolves to **v1.0, Dec 2024**, which holds
*indicative/stylized* trajectories that do NOT match the published Table 3
storylines (its 2100 CO2 for H is 74.1 GtCO2/yr vs 55.0 in the final data).
v1.1.1 carries the final IAM-based runs -- its `model` column lists the same
seven IAMs as our harmonized IIASA pull. ⚠ It does NOT simply reproduce the
prerelease: H, HL and M agree to <=0.5% on headline species, but MEDIUM-TO-LOW
(COFFEE 1.6) was substantially REVISED before publication (2100: AFOLU 55%,
Sulfur 15%, CH4 8%, N2O 9% different; total ERF up to 0.36 W/m2 / 8.5% lower
mid-century, -3.7% at 2100), and Low / Low-to-Negative / Very Low moved 1-2% in
ERF. An earlier note here claimed "<=0.3% on every headline species" -- that had
been verified for H ONLY (measured 2026-09-17). Use v1.1.1, not v1.0.

SPLICE VARIANTS (both written; they answer slightly different questions)
-----------------------------------------------------------------------
Our 1750.5-2100.5 block is IDENTICAL in both variants to the existing spliced
CSVs, so a run against the constant-tail cube is a controlled A/B on the tail.

  "harmonized"  PRIMARY. Exactly continuous at 2100 by construction, with the
             rule chosen per species by whether the species can change sign:

               CO2 FFI, CO2 AFOLU (SIGNED_SPECIES) -> ADDITIVE offset,
                 tail = ours(2100) + [ext(t) - ext(2100)].
                 These cross zero post-2100 in most markers, where a ratio is
                 undefined or sign-flipping; additive also reproduces the
                 storyline's ABSOLUTE statements directly (AFOLU ramped to zero
                 by 2150, CO2 plateau slightly above 2100, net-zero CO2 dates,
                 the "-22 GtCO2/yr by 2150" style targets in Table 3).

               every other species -> MULTIPLICATIVE ratio,
                 tail = ours(2100) * ext(t)/ext(2100).
                 These are non-negative and their extension rule is a RELAXATION
                 toward a floor ("a global relaxation level consistent with each
                 species' lowest emission value", Sect. 4), so the decay SHAPE is
                 the meaningful content and a ratio preserves it while keeping
                 the series non-negative. An additive offset here drives trace
                 F-gases and halons negative wherever the two data vintages
                 disagree on the 2100 level (up to 3662 clipped cells in VL),
                 which is an artifact of the harmonization, not of the scenario.

  "fullnative"  CALIBRATION-DECOMPOSITION variant. NOT a tail choice at all: the
             ENTIRE 1750.5-2300.5 series is taken from the extension file,
             including ITS history (Nicholls et al. 2025) rather than Smith 2024.
             Used to answer "is our ERF lower than the paper's Fig. 2 because of
             emissions, or because of the calibration vintage?" -- running our
             FaIR 1.4.5 on this reproduces OUR number, not theirs, isolating the
             gap to fair-calibrate 1.4.5-vs-1.5.0. Note this variant writes the
             complete year grid INCLUDING 2021.5 and 2022.5, which the existing
             spliced CSVs omit (SMITH_LAST_YEAR=2022 vs Smith ending 2020.5);
             writing those columns onto the old template appends them AFTER
             2300.5 and FaIR then raises NonMonotonicError.

  "zenodo"   PUBLICATION basis (2026-09-17). The ENTIRE
             scenario period 2023.5-2300.5 -- 21st century AND extension -- is
             taken verbatim from the published Zenodo v1.1.1 file, with our
             Smith 2024 history kept for 1750.5-2020.5. One published, openly
             accessible, DOI-citable source for every scenario year, and no
             harmonization at 2100 because 2023-2300 is one continuous series.
             Chosen over the IIASA prerelease (the "harmonized" variant's
             2023-2100 block) because the prerelease is GATED and is NOT the
             final record: Medium-to-Low was revised materially before
             publication (see EXTENSION SOURCE above), so a figure drawn from
             the prerelease shows a superseded ML line. H/HL/M move <0.5%.
             The splice point is the SAME one every other variant uses
             (Smith -> scenario at 2023.5), and the column set matches the
             template exactly, so the pre-2023 portion is identical to the
             other variants'.

  "native"   SENSITIVITY. Tail = ext(t) verbatim, unadjusted. Leaves a small
             step at 2100->2101 wherever the two vintages differ (<=0.3% on
             headline species; larger relative -- but negligible in forcing --
             on a few Montreal-Protocol gases whose infilling differs).

Agreement between the two is the check that the harmonization choice does not
carry the result.

UNITS: the extension file uses `Mt CO2/yr` and `kt N2O/yr`; our canonical Smith
format uses `Gt CO2/yr` and `Mt N2O/yr`. Same two conversions already applied to
the IIASA pull, reused here via EXT_TO_SMITH_SCALE.

Output: data/vanvuuren/spliced_ext_<variant>/<scenario_slug>.csv, same wide
format, same species order, and the same `scenario` column values as the
existing spliced CSVs -- so every downstream driver reads them unchanged.

PORT (Ladrillo tools/forcing/, 2026-10-08). Ported from FaIRtoFrEDI
scripts/build_vanvuuren_extended_splice.py @ 8984a13. The splice is copied
UNCHANGED; only the I/O plumbing differs:
  * the template directory (research data/vanvuuren/spliced/, built by
    build_vanvuuren_spliced.py) is --spliced-dir, and the variants are written to
    <--out-root>/spliced_ext_<variant>/ (was data/vanvuuren/).
  * --variants selects which of the four to write (default: all four, as
    before). Each variant re-reads its template, so the variants are
    independent and writing a subset changes no file.

--skeleton (Ladrillo, 2026-10-08): THE `zenodo` VARIANT WITHOUT THE TEMPLATE.
For `zenodo` every value from 2023.5 on comes from Zenodo, and the one consumer,
build_emissions_v160_cmip7harm_vv.py, reads only the 2023.5 anchor and the later
columns (its history is the CMIP7 record). The template then contributes only
its row order, the `scenario` / `region` / `variable` / `unit` strings and its
column grid -- but building it needs the gated IIASA prerelease file. So
`--skeleton vv_template_skeleton.csv` (names and units only, one row per marker
x species in the template's row order; no numbers) replaces `--spliced-dir`, and
is accepted only with `--variants zenodo`:
  * the column grid is the template's, gap included: 1750.5-2020.5, then
    2023.5-2300.5 (2021.5 and 2022.5 are ABSENT, as in the template --
    build_vanvuuren_spliced.py keeps Smith 2024 only through 2020.5);
  * EVERY column is filled from Zenodo x EXT_TO_SMITH_SCALE, so the history
    columns (1750.5-2020.5) carry Zenodo's OWN history, not Smith 2024. They
    are unused downstream: the cmip7harm rebase replaces 1750.5-2023.5 with the
    CMIP7 record and reads only 2023.5 and later from this file;
  * a species the skeleton names but Zenodo lacks is an ERROR (the template
    path keeps the template's values instead; here there are none to keep).
The columns from 2023.5 on are the same numbers the template path writes, so the
cmip7harm van Vuuren emissions are byte-identical either way (checked on all 7
markers, 2026-10-08).
Normally driven by build_forcing.py, which writes only `zenodo`, from the
skeleton -- the variant the cmip7harm van Vuuren emissions are built from.
"""

from pathlib import Path
import argparse
import sys

import numpy as np
import pandas as pd

# ── Provenance / configuration constants ─────────────────────────────────────
EXT_SOURCE_LABEL = ("van Vuuren et al. 2026 GMD 19:2627-2656 ScenarioMIP-CMIP7 "
                    "extensions; Zenodo 10.5281/zenodo.20713982 v1.1.1 "
                    "(2026-06-16) data/fair-inputs/emissions_1750-2500.csv")
OUR_SOURCE_LABEL = ("IIASA scenariomip-2026-prerelease harmonized+infilled "
                    "2023-2100, Smith 2024 history to 2022.5")

SPLICE_YEAR   = 2100      # last year taken from OUR harmonized block
TAIL_END_YEAR = 2300      # RATIFIED end year for all this work (Marcus
                          # 2026-08-28). The van Vuuren extension data runs to
                          # 2500 and is truncated here DELIBERATELY -- do not
                          # raise this because the source allows it.
VARIANTS      = ("harmonized", "native", "fullnative", "zenodo")
SCENARIO_FIRST_YEAR = 2023   # first scenario year in the template (Smith ends 2020.5)

# --skeleton: the template's column grid, reproduced. The template keeps Smith
# 2024 history 1750.5-2020.5 and starts the scenario at SCENARIO_FIRST_YEAR, so
# 2021.5 and 2022.5 are ABSENT (fair interpolates across the gap).
TEMPLATE_FIRST_YEAR = 1750
TEMPLATE_HIST_LAST_YEAR = 2020
META_COLS = ["scenario", "region", "variable", "unit"]
SKELETON_COLS = ["slug", *META_COLS]
SKELETON_VARIANT = "zenodo"  # the only variant whose values all come from Zenodo

SPLICED_DIR = None      # set from --spliced-dir in main()


def template_grid():
    """The template's year columns: 1750.5-2020.5, then 2023.5-2300.5."""
    return ([f"{y}.5" for y in range(TEMPLATE_FIRST_YEAR, TEMPLATE_HIST_LAST_YEAR + 1)]
            + [f"{y}.5" for y in range(SCENARIO_FIRST_YEAR, TAIL_END_YEAR + 1)])


def load_skeleton(path):
    """Names and units only: one row per marker x species, in template row order."""
    sk = pd.read_csv(path, dtype=str, keep_default_na=False)
    if list(sk.columns) != SKELETON_COLS:
        sys.exit(f"ERROR: {path} columns {list(sk.columns)} != {SKELETON_COLS}")
    if (sk == "").any().any():
        sys.exit(f"ERROR: {path} has empty cells")
    if set(sk["slug"]) != set(EXT_SCEN_TO_SLUG.values()):
        sys.exit(f"ERROR: {path} slugs {sorted(set(sk['slug']))} != "
                 f"{sorted(EXT_SCEN_TO_SLUG.values())}")
    if sk.duplicated(["slug", "variable"]).any():
        sys.exit(f"ERROR: {path} repeats a species within a marker")
    if (sk.groupby("slug")["scenario"].nunique() != 1).any():
        sys.exit(f"ERROR: {path} gives one marker more than one scenario name")
    return sk


def from_skeleton(rows, sub, slug):
    """One marker's wide CSV from skeleton rows + Zenodo values on the template grid."""
    cols = template_grid()
    absent = [c for c in cols if c not in sub.columns]
    if absent:
        sys.exit(f"ERROR: extension CSV lacks columns {absent[:5]}...")
    missing = [v for v in rows["variable"] if v not in sub.index]
    if missing:
        sys.exit(f"ERROR: [{slug}] species {missing} absent from Zenodo; --skeleton has "
                 f"no template values to keep -- use --spliced-dir")
    vals = np.vstack([sub.loc[spc, cols].to_numpy(dtype=float)
                      * EXT_TO_SMITH_SCALE.get(spc, 1.0) for spc in rows["variable"]])
    if np.isnan(vals).any():
        sys.exit(f"ERROR: [{slug}] NaN in the Zenodo values on the template grid")
    return pd.concat([rows[META_COLS].reset_index(drop=True),
                      pd.DataFrame(vals, columns=cols)], axis=1)

# Extension-file scenario code -> our scenario slug (= existing spliced filename)
EXT_SCEN_TO_SLUG = {
    "H":  "high_ssp3_marker",
    "HL": "high_to_low_ssp5_marker",
    "M":  "medium_ssp2_marker",
    "ML": "medium_to_low_ssp2_marker",
    "L":  "low_ssp2_marker",
    "LN": "low_to_negative_ssp2_marker",
    "VL": "very_low_ssp1_marker",
}

# Extension-file units -> canonical Smith units (see module docstring).
EXT_TO_SMITH_SCALE = {
    "CO2 FFI":   1e-3,   # Mt CO2/yr -> Gt CO2/yr
    "CO2 AFOLU": 1e-3,   # Mt CO2/yr -> Gt CO2/yr
    "N2O":       1e-3,   # kt N2O/yr -> Mt N2O/yr
}

# Species allowed to be negative (net-negative fluxes are physical for these).
# These take the ADDITIVE rule; everything else takes the MULTIPLICATIVE rule.
SIGNED_SPECIES = {"CO2 FFI", "CO2 AFOLU"}

# Ratio guard: below this the 2100 denominator is treated as zero and the
# species is held at our 2100 value (reported, never silent).
RATIO_EPS = 1e-12


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ext-csv", required=True,
                    help="path to emissions_1750-2500.csv from Zenodo 20713982 v1.1.1")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--spliced-dir", type=Path,
                     help="the constant-tail template CSVs (written by "
                          "build_vanvuuren_spliced.py; needs the gated IIASA input)")
    src.add_argument("--skeleton", type=Path,
                     help=f"names-and-units CSV ({', '.join(SKELETON_COLS)}) in place "
                          f"of the template; only with --variants {SKELETON_VARIANT}")
    ap.add_argument("--out-root", type=Path, required=True,
                    help="variants go to <out-root>/spliced_ext_<variant>/")
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS),
                    choices=VARIANTS)
    args = ap.parse_args()
    global SPLICED_DIR
    SPLICED_DIR = args.spliced_dir
    variants = [v for v in VARIANTS if v in args.variants]   # keep VARIANTS order
    skel = None
    if args.skeleton is not None:
        if variants != [SKELETON_VARIANT]:
            sys.exit(f"ERROR: --skeleton builds only --variants {SKELETON_VARIANT} "
                     f"(every other variant takes values from the template); "
                     f"got {variants}")
        skel = load_skeleton(args.skeleton)

    ext = pd.read_csv(args.ext_csv)
    if "scenario" not in ext.columns:
        sys.exit("ERROR: extension CSV has no 'scenario' column")
    ext_scens = set(ext["scenario"].unique())
    missing = set(EXT_SCEN_TO_SLUG) - ext_scens
    if missing:
        sys.exit(f"ERROR: extension CSV missing scenarios {sorted(missing)}")

    print("── van Vuuren PUBLISHED-EXTENSION splice ─────────────────────────────")
    print(f"  extension source : {EXT_SOURCE_LABEL}")
    if skel is None:
        print(f"  21st-century src : {OUR_SOURCE_LABEL}")
        print(f"  splice at        : {SPLICE_YEAR}.5   tail end: {TAIL_END_YEAR}.5")
    else:
        print(f"  names and units  : {args.skeleton} (no template; every value "
              f"from Zenodo, history columns included)")
        print(f"  grid             : {TEMPLATE_FIRST_YEAR}.5-{TEMPLATE_HIST_LAST_YEAR}.5, "
              f"{SCENARIO_FIRST_YEAR}.5-{TAIL_END_YEAR}.5 (the template's)")
    print(f"  variants         : {', '.join(variants)}")

    tail_years = list(range(SPLICE_YEAR + 1, TAIL_END_YEAR + 1))
    tail_cols  = [f"{y}.5" for y in tail_years]
    # Explicit, sorted, COMPLETE grid for the fullnative variant.
    full_years = list(range(1750, TAIL_END_YEAR + 1))
    full_cols  = [f"{y}.5" for y in full_years]

    for variant in variants:
        out_dir = args.out_root / f"spliced_ext_{variant}"
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n=== variant '{variant}' -> {out_dir} ===")

        for code, slug in EXT_SCEN_TO_SLUG.items():
            if skel is not None:
                sub = ext[ext["scenario"] == code].set_index("variable")
                built = from_skeleton(skel[skel["slug"] == slug], sub, slug)
                out = out_dir / f"{slug}.csv"
                built.to_csv(out, index=False)
                print(f"    wrote {out.name}  [{len(built)} species; every column from "
                      f"Zenodo v1.1.1; names/units from the skeleton]")
                continue

            ours = pd.read_csv(SPLICED_DIR / f"{slug}.csv")
            sub  = ext[ext["scenario"] == code].set_index("variable")

            if variant == "zenodo":
                # Every scenario-period column the template carries, from Zenodo.
                scen_cols = [c for c in ours.columns if c.endswith(".5")
                             and float(c) >= SCENARIO_FIRST_YEAR]
                n_abs = 0
                for i, spc in enumerate(ours["variable"].values):
                    if spc not in sub.index:
                        n_abs += 1
                        print(f"    [{slug}] {spc!r} absent from Zenodo -> template kept")
                        continue
                    scale = EXT_TO_SMITH_SCALE.get(spc, 1.0)
                    ours.loc[i, scen_cols] = sub.loc[spc, scen_cols].to_numpy(dtype=float) * scale
                out = out_dir / f"{slug}.csv"
                ours.to_csv(out, index=False)
                print(f"    wrote {out.name}  [{scen_cols[0]}..{scen_cols[-1]} from Zenodo v1.1.1"
                      f"{f'; {n_abs} species kept from template' if n_abs else ''}]")
                continue

            if variant == "fullnative":
                rows = []
                for _, tr in ours.iterrows():
                    spc = tr["variable"]
                    r = {k: tr[k] for k in ("scenario", "region", "variable", "unit")}
                    if spc in sub.index:
                        scale = EXT_TO_SMITH_SCALE.get(spc, 1.0)
                        vals = sub.loc[spc, full_cols].to_numpy(dtype=float) * scale
                    else:
                        print(f"    [{slug}] {spc!r} absent from extension "
                              f"-> template values retained")
                        vals = np.array([float(tr[c]) if c in ours.columns else np.nan
                                         for c in full_cols])
                    r.update(dict(zip(full_cols, vals)))
                    rows.append(r)
                out = out_dir / f"{slug}.csv"
                pd.DataFrame(rows, columns=["scenario", "region", "variable",
                                            "unit", *full_cols]).to_csv(out, index=False)
                print(f"    wrote {out.name}  [full 1750.5-{TAIL_END_YEAR}.5 from extension]")
                continue

            n_clipped = 0
            n_heldflat = 0
            max_step  = 0.0
            max_step_sp = ""
            for i, spc in enumerate(ours["variable"].values):
                if spc not in sub.index:
                    # Species absent from the extension file: keep the constant
                    # tail already present in `ours` (documented, not silent).
                    print(f"    [{slug}] {spc!r} absent from extension "
                          f"-> constant-2100 tail retained")
                    continue
                scale = EXT_TO_SMITH_SCALE.get(spc, 1.0)
                e_tail = sub.loc[spc, tail_cols].to_numpy(dtype=float) * scale
                e_2100 = float(sub.loc[spc, f"{SPLICE_YEAR}.5"]) * scale
                o_2100 = float(ours.at[i, f"{SPLICE_YEAR}.5"])

                if variant == "harmonized":
                    if spc in SIGNED_SPECIES:
                        new = o_2100 + (e_tail - e_2100)
                    elif abs(e_2100) > RATIO_EPS:
                        new = o_2100 * (e_tail / e_2100)
                    else:
                        new = np.full(len(tail_cols), o_2100)
                        n_heldflat += 1
                        print(f"    [{slug}] {spc!r} ext(2100)==0 "
                              f"-> held flat at our 2100 value")
                else:
                    new = e_tail.copy()
                    step = abs(e_2100 - o_2100)
                    if step > max_step:
                        max_step, max_step_sp = step, spc

                if spc not in SIGNED_SPECIES:
                    neg = new < 0
                    if neg.any():
                        n_clipped += int(neg.sum())
                        new = np.clip(new, 0.0, None)

                ours.loc[i, tail_cols] = new

            out = out_dir / f"{slug}.csv"
            ours.to_csv(out, index=False)
            note = f"  clipped {n_clipped} negative cells" if n_clipped else ""
            if variant == "native" and max_step_sp:
                note += f"  max 2100 step {max_step:.4g} ({max_step_sp})"
            print(f"    wrote {out.name}{note}")

    if skel is not None:
        print(f"\nDone. {SKELETON_VARIANT} from the skeleton: {SCENARIO_FIRST_YEAR}.5 "
              f"onward is the published Zenodo series; the history columns are "
              f"Zenodo's own (not Smith 2024) and are unused downstream.")
        return
    tail_only = [v for v in variants if v in ("harmonized", "native")]
    scen_from_ext = [v for v in variants if v == "zenodo"]
    print(f"\nDone. Sanity vs the template (--spliced-dir):")
    if tail_only:
        print(f"  {', '.join(tail_only)}: byte-identical through {SPLICE_YEAR}.5 "
              f"(only the post-{SPLICE_YEAR} tail differs)")
    if scen_from_ext:
        print(f"  {', '.join(scen_from_ext)}: byte-identical through "
              f"{SCENARIO_FIRST_YEAR - 1}.5 only; {SCENARIO_FIRST_YEAR}.5 onward is "
              f"the published Zenodo series")
    if "fullnative" in variants:
        print(f"  fullnative: deliberately replaces the history too -- nothing identical")


if __name__ == "__main__":
    main()
