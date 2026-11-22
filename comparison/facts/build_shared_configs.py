#!/usr/bin/env python3
"""
Write one FACTS config.yml per shared-climate experiment (10: 7 van Vuuren markers + 3 SSP controls).

Every scenario-dependent field is READ FROM THE INPUT NetCDF's OWN ATTRIBUTES, never re-typed
here: `Scenario` -> global-options.scenario and the three input filenames, `LWSScenario` ->
the landwaterstorage population pathway. A config therefore cannot disagree with the climate
file it is paired with, and a marker missing an attribute RAISES rather than defaulting
(memory: intersect_is_a_silent_default -- absence must assert, not omit).

WORKFLOWS: wf1f / wf2f / wf3f / wf4 on every key; wf1e / wf2e / wf3e (emulandice) on the ten
  BASE keys (3 SSP controls + 7 van Vuuren scenarios), 2100-capped.
  ⚠ CORRECTED 2026-09-14. Until then the "e" workflows ran on the SSP controls only, on the
  reasoning that "emulandice picks per-SSP-trained GP emulators from the scenario label". THAT
  IS NOT WHAT THE MODULE DOES: modules/emulandice/shared/emulandice_preprocess.py writes every
  forcing row with the literal scenario "FACTS", and the R side (emulandice/R/main.R, dataset
  == "FACTS") runs the GP emulators on the supplied 2015-2100 GSAT trajectories as one
  pseudo-scenario; the NetCDF Scenario attribute reaches only the output metadata. So the
  emulators take any GSAT path. The van Vuuren 2100 GSATs sit inside the ISMIP6/GlacierMIP2
  training range (SSP1-2.6 .. SSP5-8.5), so this is interpolation. The 2100 cap is real
  (emulandice years are hard-coded 2015-2100), hence no "e" workflow on any *2300 twin or
  pulse key -- their Scenario labels differ, so the set below excludes them by construction.
  ** The predecessor config (global.coupling.ssp245.fairv145.noemu) dropped the emulandice
  MODULES but left "wf1e"/"wf2e"/"wf3e" in the surviving modules' include_in_workflow. FACTS
  duly emitted total.workflow.wf1e/wf2e/wf3e files -- totals with NO ICE SHEETS IN THEM,
  indistinguishable by filename from real ones. Those strings are stripped here. **

HORIZON: pyear_end 2150 by DEFAULT, the convention every prior FACTS run on this project used
  and the horizon the shipped comparison table reports FACTS at. A climate file may carry its
  own `PyearEnd` attribute (the `*2300` twin keys: matched-dT pair 2026-09-02, seven van Vuuren
  markers 09-03, three SSP controls 09-11), which this builder reads -- so the horizon lives in
  ONE place and a config cannot disagree with its climate. The 2150 configs are never edited in
  place; the twins are separate experiments, and extract_facts_shared_components.py's [EXTEND]
  gate requires each twin to reproduce its base run bit-for-bit on the overlap before merging.
"""
import os, re, sys, glob
import numpy as np
import h5py

EXP_GLOB   = "experiments/global.shared.*.n200"
## ⚠ THE WORKFLOW SET IS PER-KEY, because emulandice is 2100-CAPPED (not because it is
## scenario-trained -- see the module docstring). The "f" workflows run everywhere; the "e"
## workflows run on the ten BASE (2150-horizon) keys and on no *2300 twin or pulse key.
## The three SSP CONTROLS -- the only keys with a prior result to check against.
CONTROLS = {"ssp126", "ssp245", "ssp585"}
VV_SCENARIOS = {"vvVL", "vvLN", "vvL", "vvML", "vvM", "vvHL", "vvH"}   # the seven van Vuuren keys
EMU_KEYS = CONTROLS | VV_SCENARIOS                 # emulandice runs here, 2100-capped
WORKFLOWS_F = ["wf1f", "wf2f", "wf3f", "wf4"]      # every key
WORKFLOWS_E = ["wf1e", "wf2e", "wf3e"]             # EMU_KEYS only (emulandice, 2100-capped)
WORKFLOWS = WORKFLOWS_F + WORKFLOWS_E               # the full vocabulary the guard allows
NSAMPS     = 200
PYEAR      = dict(start=2020, end=2150, step=10)   # DEFAULT; a file may override `end`
## ⚠ The horizon is READ FROM THE CLIMATE FILE's `PyearEnd` attribute where present, exactly as
## Scenario and LWSScenario already are, so a config cannot disagree with the climate it is paired
## with. Absent the attribute it falls back to PYEAR['end'] — every pre-2026-09-02 file.
BASEYEAR   = 2005

def wf(*names):
    """Render an include_in_workflow block, refusing any workflow not in WORKFLOWS."""
    bad = [n for n in names if n not in WORKFLOWS]
    if bad:
        raise ValueError("workflow(s) %s not in WORKFLOWS %s" % (bad, WORKFLOWS))
    return "".join('            - "%s"\n' % n for n in names)

FSET = ("wf1f", "wf2f", "wf3f")
ESET = ("wf1e", "wf2e", "wf3e")


def wfs_for(key):
    """The workflows this key may carry: emulandice on the ten base keys, none on twins/pulses."""
    return WORKFLOWS_F + (WORKFLOWS_E if key in EMU_KEYS else [])

def config_for(scen, lws_ssp, files, key, pyear_end=None):
    pyear_end = pyear_end or PYEAR['end']
    gsat, ohc, clim = files
    inputs = "".join('            - "%s"\n' % f for f in (gsat, ohc, clim))
    is_ctrl = key in EMU_KEYS   # 'ctrl' historically; now = emulandice-capable key
    ALL = tuple(wfs_for(key))
    ## emulandice blocks, emitted ONLY for the SSP controls. pyear_end 2100 is emulandice's
    ## own cap, declared here rather than discovered as a hole in the table later.
    emu = "" if not is_ctrl else """
    emuAIS:
        module_set: "emulandice"
        module: "AIS"
        include_in_workflow:
%s        options:
            pyear_end: 2100

    emuGrIS:
        module_set: "emulandice"
        module: "GrIS"
        include_in_workflow:
%s        options:
            pyear_end: 2100

    emuglaciers:
        module_set: "emulandice"
        module: "glaciers"
        include_in_workflow:
%s        options:
            pyear_end: 2100
""" % (wf("wf1e"), wf(*ESET), wf(*ESET))
    return f"""# GENERATED by build_shared_configs.py -- do not hand-edit; regenerate.
# Scenario and LWS pathway are read from {clim}'s own attributes.
global-options:
    nsamps: {NSAMPS}
    scenario: {scen}
    pyear_start: {PYEAR['start']}
    pyear_end: {pyear_end}
    pyear_step: {PYEAR['step']}
    baseyear: {BASEYEAR}
    pipeline_file: 'pipeline.global.yml'

climate_step:
    dummy:
        module_set: "facts"
        module: "dummy"
        pipeline_file: 'pipeline.yml'
        input_data_file:
{inputs}        climate_output_data:
{inputs}        generates_climate_output: true

sealevel_step:
    GrIS1f:
        module_set: "FittedISMIP"
        module: "GrIS"
        include_in_workflow:
{wf(*FSET)}
    deconto21:
        module_set: "deconto21"
        module: "AIS"
        include_in_workflow:
{wf(*(("wf3e", "wf3f") if is_ctrl else ("wf3f",)))}
    bamber19:
        module_set: "bamber19"
        module: "icesheets"
        include_in_workflow:
{wf("wf4")}
    larmip:
        module_set: "larmip"
        module: "AIS"
        include_in_workflow:
{wf(*(("wf2e", "wf2f") if is_ctrl else ("wf2f",)))}
{emu}
    ar5glaciers:
        module_set: "ipccar5"
        module: "glaciers"
        options:
            gmip: 2
        include_in_workflow:
{wf(*WORKFLOWS_F)}
    ar5AIS:
        module_set: "ipccar5"
        module: "icesheets"
        pipeline_file: "pipeline.AIS.global.yml"
        include_in_workflow:
{wf("wf1f")}
    ocean:
        module_set: "tlm"
        module: "sterodynamics"
        include_in_workflow:
{wf(*ALL)}
    lws:
        module_set: "ssp"
        module: "landwaterstorage"
        options:
            scenario: "{lws_ssp}"
            dcrate_lo: -0.4
        include_in_workflow:
{wf(*ALL)}

totaling_step:
    total:
        module_set: "facts"
        module: "total"
        loop_over_workflows: true
        loop_over_scales: true
        stages:
            - workflow
"""

REFERENCE_CFG = "experiments/global.coupling.ssp245.n200/config.yml"


def gate_against_reference():
    """The generated SSP config must reproduce the REFERENCE config's module -> workflow map.

    ⚠ THIS GATE EXISTS BECAUSE THE OBVIOUS-LOOKING VERSION WAS WRONG. Giving every shared
    module `include_in_workflow: <all workflows for this key>` reads as the safe default and is
    not: `ar5glaciers` must NOT appear in the emulandice workflows, because `emuglaciers`
    supplies glaciers there. Granting it both made wf1e/wf2e/wf3e sum TWO glacier modules and
    inflated those totals by ~11 cm (+27% at ssp245@2100). Every PER-MODULE number stayed
    correct, so nothing local looked wrong; it was caught only because nine e-cells all moved
    the same way by nearly the same amount.

    The reference is a config we did not write for this purpose and never update from here, so
    this is an EXTERNAL check rather than the script agreeing with itself
    (`gate_reads_its_own_output`).
    """
    import yaml
    ref_p = os.path.join(os.path.dirname(os.path.abspath(__file__)), REFERENCE_CFG)
    if not os.path.exists(ref_p):
        print("[GATE] ⚠ no reference config at %s -- module->workflow map UNCHECKED."
              % REFERENCE_CFG)
        return
    ref = yaml.safe_load(open(ref_p))["sealevel_step"]
    gen = yaml.safe_load(config_for("ssp245", "ssp2",
                                    ["a_gsat.nc", "a_ohc.nc", "a_climate.nc"],
                                    "ssp245"))["sealevel_step"]
    bad = []
    for mod, spec in gen.items():
        if mod not in ref:
            bad.append("%s: not in the reference config" % mod)
            continue
        a = set(spec.get("include_in_workflow", []))
        b = set(ref[mod].get("include_in_workflow", []))
        if a != b:
            bad.append("%s: generated %s, reference %s" % (mod, sorted(a), sorted(b)))
    missing = [m for m in ref if m not in gen]
    if missing:
        bad.append("in the reference but NOT generated: %s" % sorted(missing))
    if bad:
        sys.exit("[GATE] generated SSP config disagrees with %s:\n   %s"
                 % (REFERENCE_CFG, "\n   ".join(bad)))
    print("[GATE] SSP module->workflow map matches %s exactly (%d modules)\n"
          % (REFERENCE_CFG, len(gen)))


def main():
    gate_against_reference()
    exps = sorted(glob.glob(EXP_GLOB))
    if not exps:
        sys.exit("no experiment dirs matching %s" % EXP_GLOB)
    print("Writing %d configs; horizon %d-%d\n  base keys (%s) -> %s\n  twins / pulses -> %s "
          "(emulandice is 2100-capped)\n"
          % (len(exps), PYEAR['start'], PYEAR['end'], ",".join(sorted(EMU_KEYS)),
             ",".join(WORKFLOWS_F + WORKFLOWS_E), ",".join(WORKFLOWS_F)))
    for d in exps:
        clim = glob.glob(os.path.join(d, "input", "*_climate.nc"))
        if len(clim) != 1:
            sys.exit("[FAIL] %s: expected exactly 1 *_climate.nc, found %d" % (d, len(clim)))
        clim = clim[0]
        with h5py.File(clim, "r") as h:
            a = h.attrs
            for req in ("Scenario", "LWSScenario"):
                if req not in a:                       # absence ASSERTS
                    sys.exit("[FAIL] %s: missing attribute %r" % (clim, req))
            scen = a["Scenario"].decode() if isinstance(a["Scenario"], bytes) else str(a["Scenario"])
            lws  = a["LWSScenario"].decode() if isinstance(a["LWSScenario"], bytes) else str(a["LWSScenario"])
            ## Optional: absent on every pre-2026-09-02 file, which then keeps PYEAR['end'].
            pyend = a.get("PyearEnd")
            ## h5py hands back a 0-d / 1-element ndarray, not an int.
            pyend = int(np.asarray(pyend).reshape(-1)[0]) if pyend is not None else None
            grp  = [k for k in h.keys() if isinstance(h[k], h5py.Group)]
        # the group key inside the file IS the scenario string the modules will look up
        if grp != [scen]:
            sys.exit("[FAIL] %s: scenario groups %s != Scenario %r" % (clim, grp, scen))
        if lws not in ("ssp1", "ssp2", "ssp3", "ssp4", "ssp5"):
            sys.exit("[FAIL] %s: LWSScenario %r not an ssp1..ssp5 pathway" % (clim, lws))
        base = os.path.basename(clim)[: -len("_climate.nc")]
        files = ["%s_gsat.nc" % base, "%s_ohc.nc" % base, "%s_climate.nc" % base]
        for f in files:
            if not os.path.exists(os.path.join(d, "input", f)):
                sys.exit("[FAIL] %s: missing input %s" % (d, f))
        out = os.path.join(d, "config.yml")
        ## The optional emulandice block leaves a blank line behind when it is empty (the
        ## markers). Collapse runs of blank lines so a marker config is BYTE-IDENTICAL to the
        ## one its committed results were produced from -- a cosmetic diff on a config that
        ## gates a 10-experiment run is a diff someone has to re-verify for nothing.
        text = re.sub(r"\n{3,}", "\n\n", config_for(scen, lws, files, scen, pyend))
        with open(out, "w") as fh:
            fh.write(text)
        print("  %-34s scenario=%-7s lws=%-5s  wf=%s" % (os.path.basename(d), scen, lws,
              ",".join(wfs_for(scen))))
    print("\n[OK] %d configs written." % len(exps))

if __name__ == "__main__":
    main()
