#!/usr/bin/env python3
"""The four-source table behind Figs. 2, 3 and 6: Ladrillo, BRICK 2.0, MAGICC-SLR and FACTS on the seven van Vuuren
CMIP7 scenarios, per component, at 2100, 2150 and 2300.

  python python/figures/vv_model_comparison.py

Writes outputs/tables/vv_model_comparison.csv (one row per source x module x scenario x component x horizon: the
median and the 5/17/83/95 % quantiles in cm, the draw count and the band basis, plus a provenance column). Reads the
joint-arm draws of scripts/project_joint.jl (Ladrillo) and scripts/brick20_joint.jl (BRICK 2.0), their gates, and the
static comparator outputs in data/comparison/ (data/comparison/README.md says what each is).

PORT of SLR-RFF-BRICK python/vv_model_comparison.py --tag=L27 (joint, tapped arm). The table code is copied unchanged;
only the I/O differs: the draws are this repository's projections, the gates file has no CONTROL row (figlib's header
says why), and the FACTS and MAGICC-SLR outputs ship as static data (the package runs neither model). Dropped from the
research script, because no paper figure reads them: the --tag and --no-tap arms, the _width / _gmst / _pathdep side
tables and their console analyses (the Ladrillo/BRICK 2.0 width ratio, the GMST peak-vs-endpoint Spearman and the
High-to-Low minus High paired path-dependence test).

SOURCES. Ladrillo and BRICK 2.0 are each on their joint arm (posterior parameters x 841 FaIR configs, the same
2014 splice, the same 1995-2014 re-reference and the same draw-to-config permutation), so their widths are like-for-like;
gate_pairing() proves the shared permutation rather than assuming it. MAGICC-SLR runs on its own climate, computed from
the van Vuuren emissions. FACTS runs on the shared FaIR driver, n=200, and is a stack of modules per component that is
never collapsed to one number; modules that are not climate-driven at a horizon are dropped here
(figlib.FACTS_CLIMATE_DRIVEN). FACTS is relative to baseyear 2005; every other source is relative to 1995-2014.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figlib as lf  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SCRIPT = "python/figures/vv_model_comparison.py"
LADRILLO_TAG = "L27"                      # the posterior; stamped as Ladrillo's `module`, which the figures check
OUT = os.path.join(lf.TABLE_DIR, "vv_model_comparison.csv")

## Marker order is COOLEST-ENDPOINT FIRST. `family` is the CMIP6 SSP each marker's socioeconomics descend from; it is
## not an equivalence (no CMIP7 marker is a CMIP6 SSP).
MARKERS = [
    ("vvVL", "Very Low",       "SSP1"),
    ("vvLN", "Low-to-Neg",     "SSP2"),
    ("vvL",  "Low",            "SSP2"),
    ("vvML", "Medium-to-Low",  "SSP2"),
    ("vvM",  "Medium",         "SSP2"),
    ("vvHL", "High-to-Low",    "SSP5"),
    ("vvH",  "High",           "SSP3"),
]
MKEY = [m[0] for m in MARKERS]
MLABEL = {m[0]: m[1] for m in MARKERS}
MFAMILY = {m[0]: m[2] for m in MARKERS}

HORIZONS = [2100, 2150, 2300]
COMPONENTS = ["glaciers", "gis", "ais", "te", "lws", "total"]
QUANTILES = [5, 17, 50, 83, 95]

## The research table's band-basis strings, unchanged: figlib.band_is_comparable reads them, and the research table
## the paper figures were drawn from carries them. "tapped" is the Greenland tipping cell (Ladrillo's projections run
## with it on); it is checked against each draws file's own provenance below, not assumed.
BASIS_LAD = "joint (Ladrillo posterior x FaIR forcing, tapped)"
BASIS_BRK = "joint (BRICK 2.0 posterior x FaIR forcing)"
BASIS_CMP = "climate + parameter"
SRC_LAD, SRC_BRK = "Ladrillo", "BRICK 2.0"
SRC_MAG, SRC_FCT = "MAGICC-SLR", "FACTS"
## What each upstream draws file must say about itself for the basis strings above to be true.
DRAWS_MUST_SAY = {SRC_LAD: ("pairing seed 2026", "tap true"), SRC_BRK: ("pairing seed 2026",)}
## SEEDS. Every provenance clause of a draws file that names a seed is copied into this table's provenance, so the file
## Figs. 2, 3 and 6 read records every seed behind it. Each source must name these seeds (each word inside one clause),
## and its clauses must be the same in all seven draws files.
SEED_RX = re.compile(r"[^|]*seed[^|]*", re.IGNORECASE)    # also NOISE_SEED, "seeded"
SEEDS_MUST_SAY = {SRC_LAD: (("paleo", "seed"), ("pairing", "seed")),
                  SRC_BRK: (("Random.seed!", "get_model"), ("pairing", "seed"))}
DRAW_COLS = ["draw", "config", "component", "horizon", "arm", "value_cm"]

MAGICC_FILE = "magicc_nauels_components_vv.csv"
FACTS_FILE = "facts_components_shared_n200.csv"
MAGICC_CSV = os.path.join(lf.COMPARISON_DIR, MAGICC_FILE)
FACTS_CSV = os.path.join(lf.COMPARISON_DIR, FACTS_FILE)
MAGICC_MODULE = "Nauels2025"
MAGICC_N = 600          # AR6 with_slr drawnset members
FACTS_HORIZONS = [2100, 2150, 2300]


def rejected_draws(marker):
    """The joint draws the record conditioning rejected (REJECTED rows of the gates file, key draw_<k>_cfg_<config>).
    Their absence from the Ladrillo arm is by design; gate_pairing() allows exactly these and no others."""
    g = pd.read_csv(lf.require(lf.gates_csv(marker), "ladrillo", marker))
    return sorted(int(k.split("_")[1]) for k in g[g.gate == "REJECTED"].key)


def gate_ladrillo(marker):
    """figlib's projection gate, held to the research table's stricter rule: a CHECK verdict also refuses the marker
    (the research reader accepted only figlib.GATE_VERDICTS_OK)."""
    g = lf.gate_ladrillo(marker)
    if g["checks"]:
        raise SystemExit("[GATE] %s: non-passing gate row(s) (verdict CHECK): %s" % (marker, g["checks"]))
    return dict(marker=marker, control=g["control"], configs_used=g["configs"], rejected=rejected_draws(marker))


def seed_clauses(prov, path, source):
    """The provenance clauses that name a seed, after checking SEEDS_MUST_SAY."""
    cl = tuple(c.strip() for c in SEED_RX.findall(prov))
    for words in SEEDS_MUST_SAY[source]:
        if not any(all(w in c for w in words) for c in cl):
            raise SystemExit("[SEED] %s: no provenance clause names %s, so the %s seeds would go unrecorded"
                             % (os.path.relpath(path, lf.REPO), " + ".join(words), source))
    return cl


def read_joint_draws(scen, model, source):
    """One run's draws and its seed clauses, after checking the file's provenance carries what DRAWS_MUST_SAY and
    SEEDS_MUST_SAY require."""
    path = lf.require(lf.product_csv(scen, model, "draws"), model, scen)
    prov = set(pd.read_csv(path, usecols=["provenance"]).provenance)
    if len(prov) != 1:
        raise SystemExit("[PROV] %s carries %d provenance strings, expected one" % (os.path.relpath(path, lf.REPO),
                                                                                   len(prov)))
    miss = [s for s in DRAWS_MUST_SAY[source] if s not in next(iter(prov))]
    if miss:
        raise SystemExit("[PROV] %s: provenance does not say %s, so the %s band basis would be false"
                         % (os.path.relpath(path, lf.REPO), miss, source))
    return path, pd.read_csv(path, usecols=DRAW_COLS), seed_clauses(next(iter(prov)), path, source)


def load_joint(path, d, marker, source):
    """Joint-arm draws -> {(component, horizon): {q, n}} plus the config sequence."""
    d = d[d.arm == "joint"]
    if d.empty:
        raise SystemExit(f"{source} {marker}: joint arm is EMPTY in "
                         f"{os.path.relpath(path, lf.REPO)}")
    bands = {}
    for (comp, hz), g in d.groupby(["component", "horizon"]):
        q = np.percentile(g.value_cm.values, QUANTILES)
        bands[(comp, int(hz))] = dict(zip(["p05", "p17", "med", "p83", "p95"], q),
                                      n=len(g))
    s0 = d[(d.component == COMPONENTS[0]) & (d.horizon == HORIZONS[0])].sort_values("draw")
    seq = dict(zip(s0.draw.astype(int), s0.config))       # draw -> config, keyed by DRAW, not position
    return bands, seq


def load_comparator(f, source, module_col, horizons, n_default=None):
    """FACTS / MAGICC-SLR bands -> {(module, component, horizon, marker): band}.

    Both arrive as already-quantiled tidy CSVs (their per-member draws are not paired with ours), so this reads
    quantiles rather than recomputing them. A marker missing from a comparator is a FAILURE, not a skip: silently
    dropping one would leave a gap on the figure that reads as "this model has nothing to say here"."""
    if not os.path.exists(f):
        raise SystemExit(f"missing {source} arm at {os.path.relpath(f, lf.REPO)} -- see data/comparison/README.md")
    d = pd.read_csv(f)
    have = set(d.scenario.astype(str))
    gap = [m for m in MKEY if m not in have]
    if gap:
        raise SystemExit(f"[GATE] {source}: no rows for marker(s) {gap} in "
                         f"{os.path.relpath(f, lf.REPO)} -- refusing to build a four-source "
                         f"table with a hole in it.")
    out = {}
    for m in MKEY:
        for comp in COMPONENTS:
            for y in horizons:
                sub = d[(d.scenario == m) & (d.component == comp) & (d.year == y)]
                if sub.empty:
                    continue
                for r in sub.itertuples():
                    mod = str(getattr(r, module_col)) if module_col else MAGICC_MODULE
                    ## Marcus 2026-09-12: FACTS modules that are not climate-driven at this
                    ## horizon are dropped here (figlib.FACTS_CLIMATE_DRIVEN).
                    if source == SRC_FCT and not lf.facts_module_ok(mod, y):
                        continue
                    out[(mod, comp, y, m)] = dict(
                        med=float(r.med), p05=float(r.p05), p17=float(r.p17),
                        p83=float(r.p83), p95=float(r.p95),
                        n=int(getattr(r, "n", n_default) or n_default or 0))
    return out


def gate_comparator_coverage(bank, source, horizons):
    """Report which (component, horizon) cells the comparator fills, and refuse a horizon it covers only PARTIALLY.
    A source wholly absent at a horizon is a stated limitation; one present for three components and missing two is
    a broken extract. For FACTS the expected set at a horizon is the components that still have a climate-driven
    module (figlib.FACTS_CLIMATE_DRIVEN)."""
    FACTS_COMP_OF = {"larmip": "ais", "ar5AIS": "ais", "deconto21": "ais", "bamber19": None,
                     "FittedISMIP": "gis", "ar5glaciers": "glaciers", "tlm": "te", "ssp-lws": "lws"}
    lines = []
    for y in horizons:
        got = sorted({c for (_mod, c, yy, _m) in bank if yy == y})
        expect = list(COMPONENTS)
        if source == SRC_FCT and lf.FACTS_CLIMATE_DRIVEN.get(int(y)) is not None:
            keep = lf.FACTS_CLIMATE_DRIVEN[int(y)]
            expect = [c for c in COMPONENTS if c != "total"
                      and any(FACTS_COMP_OF.get(m) == c for m in keep)]
        miss = [c for c in expect if c not in got]
        if got and miss:
            raise SystemExit(f"[GATE] {source} @{y}: has {got} but is MISSING {miss}. A "
                             f"partial horizon is an incomplete extract, not a limitation.")
        scoped = [c for c in COMPONENTS if c not in expect]
        if scoped:
            print(f"[FACTS-SCOPE] @{y}: {', '.join(scoped)} carried by no climate-driven FACTS "
                  f"module -- dropped by rule, not missing.")
        lines.append((y, len(got), [c for c in COMPONENTS if c not in got]))
    return lines


def gate_pairing(seq_lad, seq_brk, marker, rejected=()):
    """LIKE-FOR-LIKE IS PROVED, NOT ASSERTED. Both drivers permute draw->config with the same pairing seed, so every
    draw present in BOTH arms must carry the same config. KEYED BY DRAW, NOT BY POSITION: the Ladrillo joint arm drops
    the draws the record conditioning rejects, and a BRICK 2.0 draw may be absent from the Ladrillo arm ONLY if
    Ladrillo's gates file lists it as REJECTED."""
    common = sorted(set(seq_lad) & set(seq_brk))
    if not common:
        raise SystemExit(f"[GATE] {marker}: the two arms share no draws to pair.")
    bad = [k for k in common if seq_lad[k] != seq_brk[k]]
    if bad:
        k = bad[0]
        raise SystemExit(
            f"[GATE] {marker}: the two models' draw->config permutations DIVERGE at draw "
            f"{k} ({seq_lad[k]} vs {seq_brk[k]}). They did not see the same forcing "
            f"subset, so no width ratio between them is like-for-like.")
    unexplained = sorted(set(seq_brk) - set(seq_lad) - set(rejected))
    if unexplained:
        raise SystemExit(f"[GATE] {marker}: BRICK 2.0 draws {unexplained[:5]} are absent from the "
                         f"Ladrillo arm and are NOT record-conditioning rejections.")
    return len(common)


def band(b):
    return f"{b['med']:7.1f} [{b['p17']:6.1f},{b['p83']:6.1f}]"


def main():
    print("Four sources on the seven van Vuuren markers: Ladrillo and BRICK 2.0 (this package, joint arms), "
          "MAGICC-SLR\n(its own climate, %d members) and FACTS (shared FaIR driver, n=200, per module), the last two "
          "from data/comparison/.\n" % MAGICC_N)

    gates, LAD, BRK, SEQ, INPUTS = {}, {}, {}, {}, []
    SEEDS = {SRC_LAD: set(), SRC_BRK: set()}
    for m in MKEY:
        gates[m] = gate_ladrillo(m)
        pl, dl, kl = read_joint_draws(m, "ladrillo", SRC_LAD)
        pb, db, kb = read_joint_draws(m, "brick20", SRC_BRK)
        SEEDS[SRC_LAD].add(kl)
        SEEDS[SRC_BRK].add(kb)
        INPUTS += [pl, pb]
        LAD[m], sl = load_joint(pl, dl, m, SRC_LAD)
        BRK[m], sb = load_joint(pb, db, m, SRC_BRK)
        SEQ[m] = gate_pairing(sl, sb, m, gates[m]["rejected"])
        print(f"[GATE] {MLABEL[m]:14s} {m:5s}  Ladrillo gates all pass "
              f"({gates[m]['control']}, {gates[m]['configs_used']} configs); "
              f"draw->config permutations agree over {SEQ[m]} draws "
              f"[Ladrillo n={LAD[m][(COMPONENTS[0], HORIZONS[0])]['n']}, "
              f"BRICK 2.0 n={BRK[m][(COMPONENTS[0], HORIZONS[0])]['n']}]")

    for src, k in SEEDS.items():
        if len(k) != 1:
            raise SystemExit("[SEED] %s: the draws files of the %d scenarios record different seeds: %s"
                             % (src, len(MKEY), sorted(k)))
    seed_text = "; ".join("%s: %s" % (src, " + ".join(next(iter(SEEDS[src])))) for src in (SRC_LAD, SRC_BRK))
    print("[SEED] the same in all %d draws files of each model: %s" % (len(MKEY), seed_text))

    MAG = load_comparator(MAGICC_CSV, SRC_MAG, None, HORIZONS, n_default=MAGICC_N)
    FCT = load_comparator(FACTS_CSV, SRC_FCT, "module", FACTS_HORIZONS)
    for src, bank, hz in ((SRC_MAG, MAG, HORIZONS), (SRC_FCT, FCT, FACTS_HORIZONS)):
        for y, ngot, _miss in gate_comparator_coverage(bank, src, hz):
            print(f"[GATE] {src:11s} @{y}: {ngot}/{len(COMPONENTS)} components present")

    rows = []
    for m in MKEY:
        # our two joint arms
        for src, bands, basis in ((SRC_LAD, LAD[m], BASIS_LAD),
                                  (SRC_BRK, BRK[m], BASIS_BRK)):
            for comp in COMPONENTS:
                for y in HORIZONS:
                    b = bands.get((comp, y))
                    if b is None:
                        continue
                    rows.append(dict(source=src,
                                     module=LADRILLO_TAG if src == SRC_LAD else "BRICK2.0",
                                     marker=m, marker_label=MLABEL[m], family=MFAMILY[m],
                                     component=comp, year=y, n_draws=b["n"],
                                     med=b["med"], p05=b["p05"], p17=b["p17"],
                                     p83=b["p83"], p95=b["p95"], band_basis=basis))
        # the two comparators; FACTS contributes one row PER MODULE
        for src, bank in ((SRC_MAG, MAG), (SRC_FCT, FCT)):
            for (mod, comp, y, mk), b in bank.items():
                if mk != m:
                    continue
                rows.append(dict(source=src, module=mod,
                                 marker=m, marker_label=MLABEL[m], family=MFAMILY[m],
                                 component=comp, year=y, n_draws=b["n"],
                                 med=b["med"], p05=b["p05"], p17=b["p17"],
                                 p83=b["p83"], p95=b["p95"], band_basis=BASIS_CMP))
    df = pd.DataFrame(rows)
    ## `scenario` is an ALIAS of `marker`, written so the figure code can key both scenario sets on one column name.
    df["scenario"] = df["marker"]
    df["provenance"] = (
        "Ladrillo %s | posterior %s | Ladrillo: %s (scripts/project_joint.jl), BRICK 2.0: %s "
        "(scripts/brick20_joint.jl), joint arm, FaIR 2.2.4 (calib 1.6.0) CMIP7-basis forcing | seeds (copied from the "
        "draws files, the same in all %d scenarios of each model): %s | MAGICC-SLR: data/comparison/%s (own climate, n=%d) | FACTS: "
        "data/comparison/%s (n=200, modules filtered by figlib.FACTS_CLIMATE_DRIVEN) | quantiles %s %% "
        "(numpy.percentile, linear) | %s; FACTS rel. baseyear 2005 | seed n/a (no RNG in this script) | commit %s"
        % (SCRIPT, LADRILLO_TAG, "outputs/projections/joint_draws_<scenario>.csv",
           "outputs/projections/brick20_draws_<scenario>.csv", len(MKEY), seed_text, MAGICC_FILE, MAGICC_N, FACTS_FILE,
           "/".join(str(q) for q in QUANTILES), lf.PROJ_BASELINE, lf.commit_stamp()))
    os.makedirs(lf.TABLE_DIR, exist_ok=True)
    df.to_csv(OUT, index=False)

    ## The four-source coverage grid: what each figure horizon can and cannot show.
    print(f"\n{'=' * 96}\nFOUR-SOURCE COVERAGE (sources present per component and horizon)"
          f"\n{'=' * 96}")
    for y in HORIZONS:
        for comp in COMPONENTS:
            have = sorted(set(df[(df.year == y) & (df.component == comp)].source))
            gap = [x for x in (SRC_LAD, SRC_BRK, SRC_MAG, SRC_FCT) if x not in have]
            print("  @%d  %-9s %-46s %s"
                  % (y, comp, ", ".join(have),
                     ("MISSING: " + ", ".join(gap)) if gap else ""))

    W = 96
    print(f"\n{'=' * W}")
    print(f"Ladrillo ({LADRILLO_TAG}, tapped) vs BRICK 2.0 — {lf.PROJ_BASELINE}, median [17-83%]")
    print("Both on their joint arm, same van Vuuren forcing: widths ARE like-for-like.")
    print("=" * W)
    for y in HORIZONS:
        print(f"\n@{y}")
        for comp in COMPONENTS:
            print(f"  --- {comp} ---")
            for m in MKEY:
                a, b = LAD[m].get((comp, y)), BRK[m].get((comp, y))
                if a is None or b is None:
                    continue
                print(f"    {MLABEL[m]:14s}  Ladrillo {band(a)}   BRICK 2.0 {band(b)}")

    print(f"\nwrote {os.path.relpath(OUT, lf.REPO)}  ({len(df)} rows)")


if __name__ == "__main__":
    main()
