#!/usr/bin/env python3
"""
Extract per-COMPONENT and per-workflow-total quantile bands from the SHARED-CLIMATE FACTS
experiments (global.shared.<key>.n200) into one tidy CSV for the four-model SLR comparison.

Sibling of extract_facts_components.py, which reads the older global.coupling.<ssp>.n200 set.
Kept separate rather than folded in, because the two sets differ in more than a path: the
shared set runs the NON-EMULANDICE workflows only, and its climate is INJECTED (there is no
`temperature.fair.temperature_gsat.nc`; the driver arrives as `shared_<key>_gsat.nc`).

Values: sea_level_change in mm rel. FACTS baseyear 2005 -> reported in cm.

⚠ BASELINE. The comparison table states MAGICC/Ladrillo/BRICK rel. 1995-2014 and FACTS rel.
  baseyear 2005, treating 2005 as ~ the 1995-2014 mean. That approximation is the STANDING
  convention here (ladrillo_model_comparison.py:27) and is NOT changed by this script -- FACTS
  output begins at pyear_start 2020 and has no 1995-2014 window to re-reference against.

⚠ MODULES, NOT A MEDIAN ACROSS MODULES. FACTS modules disagree by up to 8x, so a median over
  them summarises nothing (memory: median_needs_agreement). Every row keeps its module.

The GSAT rows carry the INJECTED driver, so the shared-climate convention can be verified from
the results themselves rather than trusted: a [DRIVER] gate below checks each extracted GSAT
against the experiment's own input file.

Usage: python3 extract_facts_shared_components.py [key ...]     (default: all 10)
"""
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

FACTS   = os.path.expanduser("~/Documents/2026/CodeProjects/facts")
OUT_CSV = os.path.expanduser(
    "~/Documents/2026/CodeProjects/SLR-RFF-BRICK/outputs/facts_components_shared_n200.csv")
EXP_TMPL = "global.shared.%s.n200"
QUANTS   = [5, 17, 50, 83, 95]

# The 10 keys: 7 van Vuuren markers + 3 SSP CONTROLS. Order = cool -> warm within each set.
# Plus the 7 van Vuuren markers extended to 2300 (2026-09-03, for the L24 deliverable Fig 4),
# same pattern as the matched-dT pair's *2300 keys in build_shared_climate_nc.py: separate
# keys, existing 2150 rows untouched.
# Plus the 3 SSP controls extended to 2300 (2026-09-11, for the L24 deliverable FIG 8, which had
# no FACTS bars at 2300). No emulandice on the twins: it is 2100-capped, so a twin is NOT in
# EMU_KEYS (the earlier 'per-SSP-trained' reason was wrong -- see build_shared_configs.py).
KEYS = ["ssp126", "ssp245", "ssp585",
        "vvVL", "vvLN", "vvL", "vvML", "vvM", "vvHL", "vvH",
        "vvVL2300", "vvLN2300", "vvL2300", "vvML2300", "vvM2300", "vvHL2300", "vvH2300",
        "ssp1262300", "ssp2452300", "ssp5852300"]
CONTROLS = {"ssp126", "ssp245", "ssp585"}
VV_SCENARIOS = {"vvVL", "vvLN", "vvL", "vvML", "vvM", "vvHL", "vvH"}
EMU_KEYS = CONTROLS | VV_SCENARIOS   # emulandice ran here (2026-09-14: the seven vv keys added)
## ⭐ THE *2300 TWINS ARE MERGED INTO THEIR BASE KEY, HERE, DETERMINISTICALLY. Until 2026-09-11 the
## merge was done BY HAND on the CSV (SLR-RFF-BRICK commit a1341eb: "vvVL2300 -> vvVL, year==2300
## only"), so re-running this extractor would have silently UNDONE it -- the CSV was not
## reproducible from its own producer. Now: each twin must reproduce its base run BIT-FOR-BIT on
## every overlapping year ([EXTEND], bound 0.0 -- same driver, same module seeds, so an identity),
## and only then are its years > REPORT_MAX_BASE appended under the base's name. The 2100/2150
## rows always come from the base experiment, so the shipped column is untouched by construction.
TWIN_OF = {k: k[: -len("2300")] for k in KEYS if k.endswith("2300")}
REPORT_MAX_BASE = 2150     # the horizon of every base experiment (build_shared_configs PYEAR)

## ⚠ THE MODULE SET IS PER-KEY, mirroring build_shared_configs.py. emulandice is 2100-CAPPED,
## so it runs on the ten BASE keys (EMU_KEYS) and on no *2300 twin. Asking for an emulandice file
## on a twin would be asking for something the config never generated, and the [FAIL] below would
## be right to fire -- so the expectation, like the config, branches on EMU_KEYS.
## (Until 2026-09-14 the vv keys were excluded on a 'per-SSP-trained' reasoning that the module
## code does not support: emulandice runs on the supplied GSAT path as one pseudo-scenario.)
WORKFLOWS_F = ["wf1f", "wf2f", "wf3f", "wf4"]      # every key
WORKFLOWS_E = ["wf1e", "wf2e", "wf3e"]             # EMU_KEYS only; emulandice is 2100-capped
SL_FILES = {
    "glaciers": [("ar5glaciers",  "ar5glaciers.ipccar5.glaciers_globalsl.nc")],
    "gis":      [("FittedISMIP",  "GrIS1f.FittedISMIP.GrIS_GIS_globalsl.nc"),
                 ("bamber19",     "bamber19.bamber19.icesheets_GIS_globalsl.nc")],
    "ais":      [("ar5AIS",       "ar5AIS.ipccar5.icesheets_AIS_globalsl.nc"),
                 ("larmip",       "larmip.larmip.AIS_globalsl.nc"),
                 ("deconto21",    "deconto21.deconto21.AIS_AIS_globalsl.nc"),
                 ("bamber19",     "bamber19.bamber19.icesheets_AIS_globalsl.nc")],
    "te":       [("tlm",          "ocean.tlm.sterodynamics_globalsl.nc")],
    "lws":      [("ssp-lws",      "lws.ssp.landwaterstorage_globalsl.nc")],
    "total":    [(wf, "total.workflow.%s.global.nc" % wf) for wf in WORKFLOWS_F],
}
## emulandice adds three MODULES and three WORKFLOW totals, on EMU_KEYS only.
SL_FILES_E = {
    "glaciers": [("emuglaciers", "emuglaciers.emulandice.glaciers_globalsl.nc")],
    "gis":      [("emuGrIS",     "emuGrIS.emulandice.GrIS_globalsl.nc")],
    "ais":      [("emuAIS",      "emuAIS.emulandice.AIS_globalsl.nc")],
    "total":    [(wf, "total.workflow.%s.global.nc" % wf) for wf in WORKFLOWS_E],
}


def files_for(key):
    """Expected (component -> [(module, fragment)]) for this key."""
    out = {c: list(v) for c, v in SL_FILES.items()}
    if key in EMU_KEYS:
        for c, v in SL_FILES_E.items():
            out.setdefault(c, []).extend(v)
    return out


# Name the payload variable explicitly. Taking list(data_vars)[0] positionally picks `lat`
# out of the injected GSAT files (which carry lat/lon alongside surface_temperature) and then
# fails somewhere else entirely -- a LENGTH/ORDER-dependent interface, the class of bug
# memory `default_off_hides_the_bug` says to guard rather than index into.
PAYLOAD = ("sea_level_change", "surface_temperature", "ocean_heat_content")

def quantile_rows(path, scenario, component, module, scale, unit_note):
    d = xr.open_dataset(path)
    cand = [v for v in PAYLOAD if v in d.data_vars]
    if not cand:
        sys.exit("[FAIL] %s: none of %s present; data_vars = %s"
                 % (os.path.basename(path), PAYLOAD, list(d.data_vars)))
    v = cand[0]
    a = d[v].squeeze()
    ycoord = next(c for c in a.dims if "year" in c.lower())
    years = np.asarray(d[ycoord].values)
    sampdim = next(dm for dm in a.dims if dm != ycoord)
    vals = a.transpose(sampdim, ycoord).values.astype(float) * scale
    rows = []
    for iy, y in enumerate(years):
        col = vals[:, iy]
        col = col[np.isfinite(col)]
        if not len(col):
            continue
        q = np.percentile(col, QUANTS)
        rows.append(dict(scenario=scenario, component=component, module=module,
                         year=int(y), med=q[2], p05=q[0], p17=q[1], p83=q[3], p95=q[4],
                         n=len(col), unit=unit_note))
    d.close()
    return rows


## Map the CONFIG's module names to the component file each one contributes to a total.
## (ar5AIS also emits EAIS/WAIS and larmip also emits PEN/SMB/WAIS; only the AIS file is summed.)
CONFIG_MODULE_FILE = {
    "ar5glaciers": "ar5glaciers.ipccar5.glaciers_globalsl.nc",
    "GrIS1f":      "GrIS1f.FittedISMIP.GrIS_GIS_globalsl.nc",
    "ar5AIS":      "ar5AIS.ipccar5.icesheets_AIS_globalsl.nc",
    "larmip":      "larmip.larmip.AIS_globalsl.nc",
    "deconto21":   "deconto21.deconto21.AIS_AIS_globalsl.nc",
    "ocean":       "ocean.tlm.sterodynamics_globalsl.nc",
    "lws":         "lws.ssp.landwaterstorage_globalsl.nc",
    "emuAIS":      "emuAIS.emulandice.AIS_globalsl.nc",
    "emuGrIS":     "emuGrIS.emulandice.GrIS_globalsl.nc",
    "emuglaciers": "emuglaciers.emulandice.glaciers_globalsl.nc",
}
## bamber19 contributes BOTH ice sheets from one module, so it needs two files.
CONFIG_MODULE_FILES_MULTI = {
    "bamber19": ["bamber19.bamber19.icesheets_AIS_globalsl.nc",
                 "bamber19.bamber19.icesheets_GIS_globalsl.nc"],
}


def gate_totals_are_sums(exp, od):
    """PER-SAMPLE identity: a workflow total must equal the sum of the modules the config
    assigns to it. This is a true IDENTITY at the sample level -- FACTS' totaling step adds the
    module arrays -- unlike `sum of medians == median of total`, which is false unless the
    components are comonotonic and must never be gated (`gate_bound_matches_its_claim`).

    ⚠ WRITTEN AFTER A DOUBLE-COUNT GOT THROUGH. Adding emulandice to the SSP controls, I also
    granted `ar5glaciers` the e-workflows; wf1e/wf2e/wf3e then summed BOTH ar5glaciers and
    emuglaciers, inflating those totals ~11 cm (+27%). Every per-module number was correct, so
    nothing local looked wrong. This gate reads the CONFIG's own composition and would have
    failed on the first run.
    """
    import yaml
    cfg_p = os.path.join(os.path.dirname(od), "config.yml")
    if not os.path.exists(cfg_p):
        return ["%s: no config.yml -- composition UNCHECKED" % exp]
    sea = yaml.safe_load(open(cfg_p))["sealevel_step"]
    wf_mods = {}
    for mod, spec in sea.items():
        for w in spec.get("include_in_workflow", []):
            wf_mods.setdefault(w, []).append(mod)
    notes = []
    for w, mods in sorted(wf_mods.items()):
        tot_p = os.path.join(od, "%s.total.workflow.%s.global.nc" % (exp, w))
        if not os.path.exists(tot_p):
            continue
        dt = xr.open_dataset(tot_p)
        tot = dt["sea_level_change"].squeeze().transpose("samples", "years").values.astype(float)
        yrs = dt["years"].values
        acc = np.zeros_like(tot)
        used = []
        for m in mods:
            frags = CONFIG_MODULE_FILES_MULTI.get(m) or [CONFIG_MODULE_FILE.get(m)]
            for fr in frags:
                if fr is None:
                    notes.append("%s %s: module %r has no file mapping" % (exp, w, m)); continue
                fp = os.path.join(od, "%s.%s" % (exp, fr))
                if not os.path.exists(fp):
                    continue
                dm = xr.open_dataset(fp)
                a = dm["sea_level_change"].squeeze().transpose("samples", "years")
                a = a.reindex(years=yrs).values.astype(float)
                acc = acc + np.nan_to_num(a)
                used.append(m)
                dm.close()
        d = np.nanmax(np.abs(np.nan_to_num(tot) - acc))
        dt.close()
        # bound: float32 storage on values up to ~1e4 mm -> a few 1e-2 mm is representation, not error
        BOUND_MM = 1.0
        status = "OK" if d <= BOUND_MM else "*** MISMATCH ***"
        if d > BOUND_MM:
            notes.append("%s %s: total - sum(%s) = %.3f mm  %s"
                         % (exp, w, "+".join(sorted(set(used))), d, status))
    return notes


def gate_twin_extends_base(twin, base):
    """[EXTEND]: the *2300 twin must equal its base experiment on every overlapping year, in
    every component file, bit-for-bit. Same injected driver, same modules, same seeds -- the
    only difference is pyear_end -- so the bound is 0.0, an identity. A non-zero here means
    the tail is NOT a continuation of the shipped 2150 run and must not be merged under its name.
    Returns a list of failure notes (empty = PASS) and the count of files compared."""
    exp_t, exp_b = EXP_TMPL % twin, EXP_TMPL % base
    od_t = os.path.join(FACTS, "experiments", exp_t, "output")
    od_b = os.path.join(FACTS, "experiments", exp_b, "output")
    notes, n = [], 0
    for comp, mods in files_for(twin).items():
        for mod, frag in mods:
            pt = os.path.join(od_t, "%s.%s" % (exp_t, frag))
            pb = os.path.join(od_b, "%s.%s" % (exp_b, frag))
            if not (os.path.exists(pt) and os.path.exists(pb)):
                notes.append("%s/%s: file missing on one side" % (twin, frag))
                continue
            dt, db = xr.open_dataset(pt), xr.open_dataset(pb)
            v = next(x for x in PAYLOAD if x in dt.data_vars)
            yc = next(c for c in dt[v].dims if "year" in c.lower())
            common = np.intersect1d(dt[yc].values, db[yc].values)
            a = dt[v].sel({yc: common}).values
            b = db[v].sel({yc: common}).values
            d = float(np.nanmax(np.abs(a - b))) if a.size else 0.0
            n += 1
            if d != 0.0 or a.shape != b.shape:
                notes.append("%s vs %s %s: max|twin-base| = %.3e over %d-%d"
                             % (twin, base, frag, d, common.min(), common.max()))
            dt.close(); db.close()
    return notes, n


def merge_twins(df, keys):
    """Append each twin's years > REPORT_MAX_BASE under its base scenario name and drop the
    twin's own rows. Only for twins whose base is ALSO in this extraction -- a twin extracted
    alone stays under its own name (there is no base column to extend)."""
    out = [df[~df.scenario.isin(TWIN_OF)]]
    n_added = 0
    for twin, base in TWIN_OF.items():
        if twin not in keys:
            continue
        tail = df[(df.scenario == twin) & (df.year > REPORT_MAX_BASE)].copy()
        if base in keys:
            tail["scenario"] = base
            n_added += len(tail)
        else:
            tail = df[df.scenario == twin]      # keep it whole, under its own name
        out.append(tail)
    return pd.concat(out, ignore_index=True), n_added


def main():
    keys = sys.argv[1:] or KEYS
    bad = [k for k in keys if k not in KEYS]
    if bad:
        sys.exit("[FAIL] unknown key(s) %s; known: %s" % (bad, KEYS))

    all_rows, driver_checks, missing, composition = [], [], [], []
    for key in keys:
        exp = EXP_TMPL % key
        od = os.path.join(FACTS, "experiments", exp, "output")
        if not os.path.isdir(od):
            missing.append("%s: no output dir" % exp)
            continue
        found = 0
        for comp, mods in files_for(key).items():
            for mod, frag in mods:
                p = os.path.join(od, "%s.%s" % (exp, frag))
                if not os.path.exists(p):
                    missing.append("%s: missing %s" % (exp, frag))
                    continue
                all_rows += quantile_rows(p, key, comp, mod, 0.1, "cm rel 2005")
                found += 1

        # the INJECTED driver, copied through by the dummy climate step
        gsat = os.path.join(od, "shared_%s_gsat.nc" % key)
        if os.path.exists(gsat):
            all_rows += quantile_rows(gsat, key, "gsat", "fair224-calib160", 1.0,
                                      "degC rel 1850-1900 (injected)")
            found += 1
            # [DRIVER] the output copy must equal the input we built
            src = os.path.join(FACTS, "experiments", exp, "input", "shared_%s_gsat.nc" % key)
            if os.path.exists(src):
                a = xr.open_dataset(gsat); b = xr.open_dataset(src)
                # same explicit payload choice as quantile_rows -- NOT list(data_vars)[0],
                # which is `lat` here and makes the gate compare NaN against NaN and pass
                # (or fail) on something that is not the climate at all.
                va = next(v for v in PAYLOAD if v in a.data_vars)
                vb = next(v for v in PAYLOAD if v in b.data_vars)
                if va != vb:
                    sys.exit("[FAIL] %s: payload differs out(%s) vs in(%s)" % (key, va, vb))
                da, db = a[va].values, b[vb].values
                if da.shape != db.shape:
                    sys.exit("[FAIL] %s: gsat shape %s != %s" % (key, da.shape, db.shape))
                if not np.isfinite(da).all():
                    sys.exit("[FAIL] %s: injected gsat has non-finite values" % key)
                dmax = float(np.max(np.abs(da - db)))
                driver_checks.append((key, dmax))
                a.close(); b.close()
        else:
            missing.append("%s: missing injected gsat" % exp)
        comp_notes = gate_totals_are_sums(exp, od)
        composition += comp_notes
        print("%-26s extracted %2d files%s"
              % (exp, found, "   [COMPOSITION FAIL]" if comp_notes else ""))

    if missing:
        print("\n[FAIL] %d missing input(s) -- refusing to write a partial table:" % len(missing))
        for m in missing:
            print("   ", m)
        sys.exit(1)

    # [COMPOSITION] gate: each workflow total == the per-sample sum of its config's modules.
    print("\n[COMPOSITION] workflow totals vs the per-sample sum of their config modules")
    if composition:
        for c in composition:
            print("   " + c)
        sys.exit("[FAIL] [COMPOSITION] %d workflow total(s) do not equal their module sum."
                 % len(composition))
    print("   -> PASS: every workflow total equals the sum of the modules its config assigns.")

    # [DRIVER] gate: the climate FACTS used must be bit-equal to the climate we handed it.
    print("\n[DRIVER] injected GSAT, output-vs-input max|diff| (bound 0.0, an IDENTITY: it is a copy)")
    worst = max(d for _, d in driver_checks) if driver_checks else 0.0
    for k, d in driver_checks:
        print("   %-8s %.3e %s" % (k, d, "OK" if d == 0.0 else "*** DIFFERS ***"))
    if worst != 0.0:
        sys.exit("[FAIL] [DRIVER] the climate FACTS ran on is not the climate we built (%.3e)" % worst)
    print("   -> PASS on all %d: FACTS ran on exactly the shared driver." % len(driver_checks))

    # [EXTEND] gate: every *2300 twin reproduces its base run on the overlap, bit-for-bit.
    twins = [k for k in keys if k in TWIN_OF and TWIN_OF[k] in keys]
    if twins:
        print("\n[EXTEND] *2300 twins vs their base experiments on the overlapping years "
              "(bound 0.0, an IDENTITY: same driver, same modules, same seeds)")
        bad = []
        for t in twins:
            notes, n = gate_twin_extends_base(t, TWIN_OF[t])
            print("   %-12s -> %-8s %2d files  %s" % (t, TWIN_OF[t], n, "OK" if not notes
                                                      else "*** DIFFERS ***"))
            bad += notes
        if bad:
            for b in bad:
                print("   " + b)
            sys.exit("[FAIL] [EXTEND] %d file(s) where the 2300 tail is not a continuation "
                     "of the shipped run; refusing to merge." % len(bad))
        print("   -> PASS on all %d: each tail continues its base run; years > %d are merged "
              "under the base name." % (len(twins), REPORT_MAX_BASE))

    df = pd.DataFrame(all_rows)
    df, n_merged = merge_twins(df, keys)
    print("[MERGE] %d twin rows (years > %d) appended under their base scenario names; "
          "twin keys dropped from the table." % (n_merged, REPORT_MAX_BASE))
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print("\nwrote %s  (%d rows, %d scenarios)" % (OUT_CSV, len(df), df.scenario.nunique()))

    for key in keys:
        for y in (2100, 2150, 2300):
            sub = df[(df.scenario == key) & (df.year == y) & (df.component != "gsat")]
            if not len(sub):
                continue
            tag = "CONTROL" if key in CONTROLS else "van Vuuren marker"
            print("\n=== %s @%d  [%s] (cm rel 2005, median [17-83%%]) ===" % (key, y, tag))
            for comp in ["glaciers", "gis", "ais", "te", "lws", "total"]:
                for _, r in sub[sub.component == comp].iterrows():
                    print("  %-9s %-12s %7.1f  [%6.1f, %6.1f]"
                          % (comp, r["module"], r["med"], r["p17"], r["p83"]))


if __name__ == "__main__":
    main()
