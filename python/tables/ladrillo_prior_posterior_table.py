#!/usr/bin/env python3
"""
ladrillo_prior_posterior_table.py -- Table A1 for the GMD paper: the prior / posterior table for every sampled
Ladrillo parameter (modelled on Wong et al. 2017 GMD Tables A1-A7).

  python python/tables/ladrillo_prior_posterior_table.py [--outdir=outputs/tables]

PORT of SLR-RFF-BRICK python/ladrillo_prior_posterior_table.py --tag=L27. The table code is copied unchanged; only
the I/O and the glacier-block names differ (the package's SLOWG / FASTG are the research code's SLOWP / FAST; R19 is
unchanged), so the research script's final label mapping (`SLOWP` -> `SLOWG`, `_FAST` -> `_FASTG`) is dropped:
the names arrive already renamed, and mapping them again would print `_FASTGG`.

Inputs
  data/constants/ladrillo_priors_L27.csv    -- the prior AS THE CALIBRATOR SCORES IT (calibrate_mcmc_ext.jl
                                               --dump-priors with the production flags; FREE, PRIOR_SKIP, the paleo
                                               block, the noise block), never transcribed by hand; imported from the
                                               research repository with the block names renamed (data/MANIFEST.csv)
  data/posterior/ladrillo_L27_posterior.csv -- the shipped posterior (10,000 draws)
Outputs
  <outdir>/ladrillo_prior_posterior.csv  (+ provenance column; research: outputs/ladrillo_prior_posterior_L27.csv)
  <outdir>/ladrillo_prior_posterior.md   (the appendix table, grouped by block; research: ..._L27.md)

DISPLAY BASE (Marcus 2026-11-10; research commit b8d7fdf): parameters the calibrator samples as log10 (LOG10_SAMPLED)
are SHOWN in the markdown table as natural logs, so that every log-scale row of Table A1 uses one base (gis_slow_ell and
ais_precip_u are ln in the code). The prior and the posterior are unchanged: the md multiplies the prior mu, sigma,
bounds and the posterior quantiles by ln 10. The CSV keeps the code's log10 values and units.

GATES (hard-fail): the prior file's parameter order == the posterior's columns (so the counts agree); every posterior
draw lies inside its prior bounds (rho < 0.99); every parameter appears in exactly one block. [LN-DISPLAY] each ln
row's quantiles, recomputed independently as quantiles of ln(10**x), equal the scaled log10 quantiles to 1e-9, and its
units string starts with "ln(". A prior row without a description entry below also fails -- a parameter must not
appear in the paper without its units.
python/tables/check_tables_identity.py gates the outputs against the research files.

Only the DESCRIPTION/UNITS column is hand-written here. Sources: DAIS units from Wong et al. 2017
Table A4; glacier-block units from julia/glaciers_nu3_component.jl (m SLE, K glacier-frame,
log10 kappa); Greenland from julia/greenland_ab_component.jl (m SLE); ledger/scope terms from
calibrate_mcmc_ext.jl comments (mm); d2 in cm (unit-RMS basis); noise in cm (residuals in cm).
(Paths are the research repository's.)
"""
import argparse
import os
import subprocess
from datetime import datetime

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = "python/tables/ladrillo_prior_posterior_table.py"
TAG = "L27"                                       # the posterior's calibration
PRIOR_FILE = os.path.join(REPO, "data", "constants", f"ladrillo_priors_{TAG}.csv")
POST_FILE = os.path.join(REPO, "data", "posterior", f"ladrillo_{TAG}_posterior.csv")
OUT_STEM = "ladrillo_prior_posterior"
DEFAULT_OUTDIR = os.path.join(REPO, "outputs", "tables")
LN10 = np.log(10.0)
# sampled as log10 by the calibrator; the md table shows them as ln (see DISPLAY BASE above). ais_ramp_log10s is not
# sampled in L27 and is kept so the set matches the research script's.
LOG10_SAMPLED = {"gic_log10_kappa_R19", "gic_log10_kappa_SLOWG", "gic_log10_kappa_FASTG", "gis_s_high",
                 "ais_ramp_log10s"}
LN_NOTE = ("ᶜ sampled as log₁₀ in the code and the posterior file; shown here as ln = ln 10 × log₁₀ "
           "(ln 10 = %.3f)." % LN10)

# block -> (heading, [parameter names in table order])
# The first heading is the paper's text, kept verbatim: its "in the code" refers to the research code's names.
BLOCKS = [
    ("Glacier reservoirs (three blocks; SLOWG = SLOWP and FASTG = FAST in the code)", [
        "gic_a_R19", "gic_b_R19", "gic_T_off_R19", "gic_log10_kappa_R19",
        "gic_a_SLOWG", "gic_b_SLOWG", "gic_T_off_SLOWG", "gic_log10_kappa_SLOWG",
        "gic_a_FASTG", "gic_b_FASTG", "gic_T_off_FASTG", "gic_log10_kappa_FASTG"]),
    ("Glacier regional temperature amplification", ["gic_amp_R19", "gic_amp_SLOWG", "gic_amp_FASTG"]),
    ("Glacier scope and ledger terms (hindcast-target constructs)", ["gic_u_unch", "gic_delta", "gic_u_pre", "gic_s_r5"]),
    ("Greenland ice sheet (two channels, two basins)", [
        "gis_c1", "gis_c0", "gis_f", "gis_alpha_f", "gis_beta_f", "gis_slow_ell", "gis_slow_w",
        "gis_s_high", "gis_amp"]),
    ("Antarctic ice sheet (DAIS) and Antarctic ocean", [
        "ais_gmst_amp", "ais_ocean_temperature₀", "antarctic_alpha", "antarctic_nu",
        "antarctic_temp_threshold", "anto_alpha", "anto_beta", "antarctic_lambda",
        "antarctic_gamma", "antarctic_kappa",
        "ais_mu", "ais_bedheight0", "ais_slope", "ais_iceflow0", "ais_precip0_LOG", "ais_precip_u",
        "ais_runoff_Ton", "ais_c",
        # L30 (--ais-ramp): present only in a ramp posterior; the BLOCKS gate below ignores names the
        # prior file does not carry, so this line is inert for L27.
        "ais_ramp_gon", "ais_ramp_log10s"]),
    ("Thermal expansion", ["thermal_alpha"]),
    ("Model-discrepancy coefficients", ["d2_gsic_1", "d2_gsic_2", "d2_steric_1", "d2_steric_2"]),
    ("Observation-error model (AR(1) per target series)", [
        "sd_ais", "rho_ais", "sd_gsic", "rho_gsic", "sd_gis", "rho_gis", "sd_steric", "rho_steric"]),
]

DESC = {
    "gic_a_R19": ("Committed-loss scale a, RGI 19", "m SLE"),
    "gic_b_R19": ("Equilibrium temperature sensitivity b, RGI 19", "K⁻¹ (glacier frame)"),
    "gic_T_off_R19": ("Equilibrium temperature offset T_off, RGI 19", "K rel. 1850–1900, glacier frame"),
    "gic_log10_kappa_R19": ("log₁₀ of the response-rate constant κ, RGI 19", "log₁₀(yr⁻¹ K⁻ν)"),
    "gic_a_SLOWG": ("Committed-loss scale a, SLOWG", "m SLE"),
    "gic_b_SLOWG": ("Equilibrium temperature sensitivity b, SLOWG", "K⁻¹ (glacier frame)"),
    "gic_T_off_SLOWG": ("Equilibrium temperature offset T_off, SLOWG", "K rel. 1850–1900, glacier frame"),
    "gic_log10_kappa_SLOWG": ("log₁₀ of the response-rate constant κ, SLOWG", "log₁₀(yr⁻¹ K⁻ν)"),
    "gic_a_FASTG": ("Committed-loss scale a, FASTG", "m SLE"),
    "gic_b_FASTG": ("Equilibrium temperature sensitivity b, FASTG", "K⁻¹ (glacier frame)"),
    "gic_T_off_FASTG": ("Equilibrium temperature offset T_off, FASTG", "K rel. 1850–1900, glacier frame"),
    "gic_log10_kappa_FASTG": ("log₁₀ of the response-rate constant κ, FASTG", "log₁₀(yr⁻¹ K⁻ν)"),
    "gic_amp_R19": ("Regional/global warming ratio, RGI 19", "–"),
    "gic_amp_SLOWG": ("Regional/global warming ratio, SLOWG", "–"),
    "gic_amp_FASTG": ("Regional/global warming ratio, FASTG", "–"),
    "gic_u_unch": ("Uncharted-ice content of the Frederikse glacier target", "mm SLE"),
    "gic_delta": ("Early-segment (1900–1960) rate bias of the glacier target", "mm yr⁻¹"),
    "gic_u_pre": ("Pre-1901 uncharted-ice set-aside (ledger)", "mm SLE"),
    "gic_s_r5": ("RGI 05 share of the 19th-century glacier datum (ledger)", "mm SLE"),
    "gis_c1": ("Committed-loss sensitivity to regional temperature", "m SLE K⁻¹"),
    "gis_c0": ("Committed loss at zero regional anomaly", "m SLE"),
    "gis_f": ("Surface-mass-balance share of the committed loss", "–"),
    "gis_alpha_f": ("Fast (SMB) channel rate, temperature-dependent part", "yr⁻¹ K⁻¹"),
    "gis_beta_f": ("Fast (SMB) channel rate, constant part", "yr⁻¹"),
    "gis_slow_ell": ("Slow (discharge) channel: log rate at the reference temperature", "ln(yr⁻¹)"),
    "gis_slow_w": ("Slow (discharge) channel: temperature-dependent fraction of the rate", "–"),
    "gis_s_high": ("log₁₀ rate scale of the high basin (NO+NE) relative to the active basin", "log₁₀(–)"),
    "gis_amp": ("Southern-Greenland/global warming ratio", "–"),
    "ais_gmst_amp": ("Antarctic/global warming ratio", "–"),
    # L30 (--ais-ramp): the additional discharge response. The onset is SAMPLED in global warming and the
    # threshold derived per draw through that draw's amp, so the table reports the sampled coordinate.
    "ais_ramp_gon": ("Onset of the additional discharge response, in global warming", "°C"),
    "ais_ramp_log10s": ("log10 slope of the additional discharge response, per °C of Antarctic excess", "log10(m SLE yr⁻¹ °C⁻¹)"),
    "ais_ocean_temperature₀": ("Initial high-latitude ocean subsurface temperature T_oc,0", "°C"),
    "antarctic_alpha": ("Partition of ocean subsurface temperature into ice flux, α_DAIS", "–"),
    "antarctic_nu": ("Runoff-decrease-with-height to precipitation constant ν", "m⁻¹ᐟ² yr⁻¹ᐟ²"),
    "antarctic_temp_threshold": ("Fast-dynamics temperature threshold", "°C (DAIS scale)"),
    "anto_alpha": ("Antarctic ocean temperature sensitivity a_ANTO", "°C °C⁻¹"),
    "anto_beta": ("Antarctic ocean temperature offset b_ANTO", "°C"),
    "antarctic_lambda": ("Fast-dynamics disintegration rate λ", "m yr⁻¹"),
    "antarctic_gamma": ("Power of ice-flow speed on water depth γ", "–"),
    "antarctic_kappa": ("Exponential dependence of precipitation on Antarctic temperature κ_DAIS", "°C⁻¹"),
    "ais_mu": ("Profile parameter µ", "m¹ᐟ²"),
    "ais_bedheight0": ("Bed height at the continent centre b₀", "m"),
    "ais_slope": ("Bed slope", "–"),
    "ais_iceflow0": ("Ice-flow constant f₀", "m yr⁻¹"),
    "ais_precip0_LOG": ("ln of the precipitation constant P₀", "ln(m yr⁻¹)"),
    "ais_precip_u": ("u = ln P₀ + κ_DAIS·T̄ (T̄ = −17.99 °C, DAIS scale); ln P₀ is derived", "ln(m yr⁻¹)"),
    "ais_runoff_Ton": ("Runoff onset temperature T_on = −h₀/c (sampled in place of h₀)", "°C (DAIS scale)"),
    "ais_c": ("Runoff-line slope c", "m °C⁻¹"),
    "thermal_alpha": ("Thermal expansion coefficient α", "kg m⁻³ °C⁻¹"),
    "d2_gsic_1": ("Glacier discrepancy, basis coefficient 1", "cm (RMS)"),
    "d2_gsic_2": ("Glacier discrepancy, basis coefficient 2", "cm (RMS)"),
    "d2_steric_1": ("Thermal-expansion discrepancy, basis coefficient 1", "cm (RMS)"),
    "d2_steric_2": ("Thermal-expansion discrepancy, basis coefficient 2", "cm (RMS)"),
    "sd_ais": ("AR(1) innovation sd, Antarctica", "cm"),
    "rho_ais": ("AR(1) autocorrelation, Antarctica", "–"),
    "sd_gsic": ("AR(1) innovation sd, glaciers", "cm"),
    "rho_gsic": ("AR(1) autocorrelation, glaciers", "–"),
    "sd_gis": ("AR(1) innovation sd, Greenland", "cm"),
    "rho_gis": ("AR(1) autocorrelation, Greenland", "–"),
    "sd_steric": ("AR(1) innovation sd, thermal expansion", "cm"),
    "rho_steric": ("AR(1) autocorrelation, thermal expansion", "–"),
}


# md-only description/units for the LOG10_SAMPLED rows (the CSV keeps DESC, which matches its log10 values)
DESC_LN = {
    "gic_log10_kappa_R19": ("ln of the response-rate constant κ, RGI 19", "ln(yr⁻¹ K⁻ν)ᶜ"),
    "gic_log10_kappa_SLOWG": ("ln of the response-rate constant κ, SLOWG", "ln(yr⁻¹ K⁻ν)ᶜ"),
    "gic_log10_kappa_FASTG": ("ln of the response-rate constant κ, FASTG", "ln(yr⁻¹ K⁻ν)ᶜ"),
    "gis_s_high": ("ln rate scale of the high basin (NO+NE) relative to the active basin", "ln(–)ᶜ"),
    "ais_ramp_log10s": ("ln slope of the additional discharge response, per °C of Antarctic excess",
                        "ln(m SLE yr⁻¹ °C⁻¹)ᶜ"),
}
assert set(DESC_LN) == LOG10_SAMPLED, "DESC_LN must cover LOG10_SAMPLED exactly"


def fmt(v, ref):
    """Format v to a precision set by ref (the parameter's scale)."""
    if v is None or (isinstance(v, float) and (np.isnan(v) or np.isinf(v))):
        return "–" if v is None or np.isnan(v) else ("∞" if v > 0 else "−∞")
    s = abs(ref) if ref not in (None, 0) and np.isfinite(ref) else abs(v)
    if s == 0:
        return "0"
    d = max(0, 2 - int(np.floor(np.log10(s))))        # three significant figures on the row's scale
    d = min(d, 6)
    t = f"{v:.{d}f}".rstrip("0").rstrip(".") if d > 0 else f"{v:.0f}"
    return t.replace("-", "−")


def prior_text(r):
    form = r.prior_form
    if form.startswith("flat"):
        return f"flat [{fmt(r.lo, r.lo)}, {fmt(r.hi, r.hi)}]"
    if form.startswith("half-normal"):
        return "half-normal N⁺(0, 5)"
    if form.startswith("N(k10c"):
        return f"N(κ̂(amp), {fmt(r.sigma, r.sigma)}) on [{fmt(r.lo, r.lo)}, {fmt(r.hi, r.hi)}]ᵃ"
    if form.startswith("joint paleo"):
        return f"paleo N({fmt(r.mu, r.mu)}, {fmt(r.sigma, r.sigma)}) on [{fmt(r.lo, r.lo)}, {fmt(r.hi, r.hi)}]ᵇ"
    return f"N({fmt(r.mu, r.mu)}, {fmt(r.sigma, r.sigma)}) on [{fmt(r.lo, r.lo)}, {fmt(r.hi, r.hi)}]"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    a = ap.parse_args()
    outdir = os.path.abspath(a.outdir)
    pri = pd.read_csv(PRIOR_FILE)
    post_path = POST_FILE
    post = pd.read_csv(post_path)
    names = list(pri.name)
    # ---- gates ----
    npost = [c for c in post.columns if c != "log_post"]
    assert len(names) == len(npost), f"prior file has {len(names)} parameters, posterior has {len(npost)}"
    assert list(npost) == names, "prior file order != posterior column order"
    missing = [n for n in names if n not in DESC]
    assert not missing, f"no description/units for {missing}"
    blocks = [(h, [n for n in ps if n in names]) for h, ps in BLOCKS]     # tag-aware: L26 drops 3, renames 1
    listed = [n for _, ps in blocks for n in ps]
    # A name listed in BLOCKS but absent from this tag's prior file is dropped (the L30 ramp pair on a
    # pre-L30 tag); a name in the prior file but NOT in BLOCKS still fails, which is the direction that matters.
    listed = [n for n in listed if n in names]
    assert sorted(listed) == sorted(names), f"BLOCKS does not cover every parameter exactly once: {sorted(set(names)-set(listed))}"
    for r in pri.itertuples():
        v = post[r.name].to_numpy()
        lo, hi = r.lo, r.hi
        if r.name.startswith("rho_"):
            hi = 0.99
        assert (v >= lo - 1e-12).all() and (v <= hi + 1e-12).all(), \
            f"{r.name}: posterior draws outside the prior bounds [{lo}, {hi}] ({v.min()}, {v.max()})"
    # ---- table ----
    q = post[names].quantile([0.05, 0.5, 0.95]).T
    q.columns = ["p05", "p50", "p95"]
    out = pri.set_index("name").join(q)
    out["mean"] = post[names].mean()
    out["description"] = [DESC[n][0] for n in out.index]
    out["units"] = [DESC[n][1] for n in out.index]
    commit = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    out["provenance"] = (f"Ladrillo {SCRIPT} | tag {TAG} | priors from calibrate_mcmc_ext.jl --dump-priors "
                         f"({pri.provenance.iloc[0].split(' | ')[1]}; {os.path.relpath(PRIOR_FILE, REPO)}) | "
                         f"posterior {os.path.relpath(post_path, REPO)} n={len(post)} | "
                         f"quantiles 5/50/95 over all draws | commit {commit} | {datetime.now():%Y-%m-%d %H:%M}")
    os.makedirs(outdir, exist_ok=True)
    csv_out = os.path.join(outdir, f"{OUT_STEM}.csv")
    out.reset_index().to_csv(csv_out, index=False)
    # ---- markdown, grouped by block ----
    lines = [f"**Table A1.** Prior distributions and posterior median (5–95%) of the {len(names)} sampled Ladrillo "
             f"parameters ({len(post):,} thinned draws of the four chains). N(µ, σ) on [lo, hi] is a normal density "
             "restricted to the bounds; 'flat' is bounds only. Glacier-block temperatures are in the block's regional "
             "frame (K relative to 1850–1900). ᵃ κ's prior centre κ̂(amp) is a log-linear function of the block's sampled "
             "amplification (τ₅₀ anchored to GlacierMIP3); ᵇ the seven DAIS geometry parameters carry a joint normal prior "
             "with the correlation of the DAIS paleo ensemble — the marginal is shown. The Greenland channel ordering "
             "(SMB faster than discharge) is a hard constraint on top of these marginals. Antarctic temperatures on the "
             "DAIS scale; ais_runoff_Ton replaces h₀ (h₀ = −T_on·c).", ""]
    lines += ["| Parameter | Description | Units | Prior | Median | 5% | 95% |",
              "|---|---|---|---|---|---|---|"]
    shown = []
    for head, ps in blocks:
        if not ps:
            continue
        lines.append(f"| **{head}** | | | | | | |")
        for n in ps:
            r = out.loc[n].copy()
            pr = pri.set_index('name').loc[n].copy()
            if n in LOG10_SAMPLED:                           # display in ln; the CSV keeps log10
                for k in ("p05", "p50", "p95"):
                    r[k] = r[k] * LN10
                for k in ("mu", "sigma", "lo", "hi"):
                    pr[k] = pr[k] * LN10
                r["description"], r["units"] = DESC_LN[n]
                lnq = np.quantile(np.log(10.0 ** post[n].to_numpy()), [0.05, 0.5, 0.95])
                assert np.allclose(lnq, [r.p05, r.p50, r.p95], rtol=0, atol=1e-9), \
                    f"[LN-DISPLAY] {n}: scaled quantiles {[r.p05, r.p50, r.p95]} != ln(10**x) quantiles {lnq}"
                assert r.units.startswith("ln("), f"[LN-DISPLAY] {n}: units '{r.units}' do not say ln"
                shown.append(n)
            code = n                                          # already the package name (SLOWG / FASTG)
            sc = max(abs(r.p05), abs(r.p50), abs(r.p95))     # one precision per row
            lines.append(f"| `{code}` | {r.description} | {r.units} | {prior_text(pr)} | "
                         f"{fmt(r.p50, sc)} | {fmt(r.p05, sc)} | {fmt(r.p95, sc)} |")
    if shown:
        lines[0] = lines[0] + " " + LN_NOTE
    md_out = os.path.join(outdir, f"{OUT_STEM}.md")
    with open(md_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n\n" + f"<!-- {out.provenance.iloc[0]} -->\n")
    print(f"wrote {csv_out} and {md_out} ({len(names)} parameters; all gates PASS)")
    print(f"[LN-DISPLAY] PASS: {len(shown)} log10-sampled rows shown as ln ({', '.join(shown)})")


if __name__ == "__main__":
    main()
