# Likelihood weights for the draws whose Antarctic fast dynamics switches on inside the observed record.
#
#   julia --project=. scripts/record_likelihood_weights.jl [--outdir=outputs/conditioning]
#
# THE WEIGHT. L27 was sampled with lambda / T_crit at their paleo medians, which never fire over the record, and the
# propagated (lambda, T_crit) pair is attached afterwards from the paleo prior. The posterior conditioned on the record is
# then reached exactly by importance weighting:
#     w_i = L(theta_i, lambda_i, Tcrit_i) / L(theta_i, lambda_med, Tcrit_med)
# On the calibration climate, a draw whose pair does not fire within the record has a hindcast identical to the
# calibration's, so w = 1 EXACTLY (asserted on all 10,000 draws, not assumed). Only the Antarctic level term can move:
# fast dynamics enters nothing but ais_sea_level, and the SMB anchor is a 1979-2008 mean, before any onset. So
#     log w = logL_ais(fast dynamics on) - logL_ais(T_crit = +Inf).
#
# THE LIKELIHOOD is the calibrator's Antarctic term (SLR-RFF-BRICK julia/calibrate_mcmc_ext.jl hetero_logl_ar1, with
# L27's --obs-corr-len=100 and the default --ais-fit-from=1900): AR(1) model error with the draw's own sd_ais and rho_ais,
# plus the per-year observational band as a correlated error with e-folding length OBS_CORR_LEN_YR,
#     Sigma = sd^2/(1-rho^2) rho^|i-j| + eps_i eps_j exp(-|i-j|/L),   eps = max(90% half-band / 1.645, 0.05 cm),
# on the residual (model - target), re-referenced to 1995-2005, over every year with an Antarctic target.
# Hindcast conventions are the calibrator's: ssp245harm mean forcing, 1850-2026, land water :central, tap off.
# [RESID] the fast-dynamics-ON residuals must equal outputs/hindcast/ic_hindcast_residuals_ladrillo.csv (the Table 5
# input, gate-4 identical to the research code) for every draw and year.
#
# JOINT ARM (a choice, flagged): the calibration likelihood is defined on the mean climate. For a joint-arm pair crossing
# on its own FaIR config the script reports the SAME ratio evaluated on that config's spliced climate: it isolates the
# fast-dynamics channel and is not an importance weight for the climate itself.
#
# NOISE-MODEL SENSITIVITY: three more lenient error models are scored alongside (dlogl_<variant>; see NOISE_VARIANTS).
# They are reported, never used for the weight.
#
# Writes <outdir>/record_likelihood_weights.csv (+ provenance column).

using Ladrillo, CSV, DataFrames, LinearAlgebra, Printf, Random, Statistics

arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const OUTDIR = arg("--outdir=", joinpath(dirname(@__DIR__), "outputs", "conditioning"))
const OBS_CORR_LEN_YR = 100.0              # L27's --obs-corr-len=100 (run_mcmc_L27.sh)
const AIS_FIT_FROM = 1900                  # the calibrator's default --ais-fit-from
const Y0, Y1 = 1850, Ladrillo.HINDCAST_END
const FORCING = Ladrillo.RECORD_FORCING    # ssp245harm, the calibration forcing
const RESID_FILE = joinpath(dirname(@__DIR__), "outputs", "hindcast", "ic_hindcast_residuals_ladrillo.csv")
const PROJ_ROWS = 1:5:10000
const PAIR_SEED = 2026                     # scripts/project_joint.jl
const RE = record_end()
const PKG_COMMIT = try readchomp(`git -C $(dirname(@__DIR__)) rev-parse --short HEAD`) catch; "no git" end

## ---- the Antarctic target, exactly as the calibrator builds it ------------------------------------------------------
tg = Ladrillo.calibration_targets(; years = AIS_FIT_FROM:Y1)
const AY = [y for (k, y) in enumerate(tg.years) if !isnan(tg.obs[:ais][k])]
AY == collect(first(AY):last(AY)) || error("Antarctic target has a year gap (AR(1) needs unit spacing)")
last(AY) == RE || error("Antarctic target ends $(last(AY)), record_end() says $RE")
const AOBS = [tg.obs[:ais][findfirst(==(y), tg.years)] for y in AY]
const AEPS = [tg.sigma[:ais][findfirst(==(y), tg.years)] for y in AY]
const HDIST = abs.(collect(1:length(AY))' .- collect(1:length(AY)))
## NOISE-MODEL SENSITIVITY (reported, not used for the weight): the band collapses from ~0.09 cm (2014-18) to the
## 0.05 cm floor at the 2019 GRACE-FO splice, and a correlated band (L = 100) cannot absorb a one-year jump, so a draw
## whose excess appears only in the last year is the case a more lenient error model could keep.
const BAND_PRE = mean(AEPS[(AY .>= 2014) .& (AY .<= 2018)])
const NOISE_VARIANTS = [
    ("l27",       AEPS, OBS_CORR_LEN_YR),                                      # the calibrator's: THE weight
    ("diag",      AEPS, 0.0),                                                  # band independent per year (L24 form)
    ("diag_wide", [y >= 2019 ? max(e, BAND_PRE) : e for (e, y) in zip(AEPS, AY)], 0.0),   # undo the 2019 collapse
    ("diag_x3",   3 .* AEPS, 0.0)]                                             # band tripled
const BANDCOVS = Dict(n => (L > 0 ? (e * e') .* exp.(-HDIST ./ L) : Matrix(Diagonal(e .^ 2))) for (n, e, L) in NOISE_VARIANTS)

"""The calibrator's hetero_logl_ar1 (variant "l27": L > 0): logpdf(MvNormal(Symmetric(Sigma)), res), written out
through a Cholesky factor (Distributions is not a dependency of this package)."""
function logl_ais(res, sd, rho, variant="l27")
    Σ = (sd^2 / (1 - rho^2)) .* rho .^ HDIST .+ BANDCOVS[variant]
    C = cholesky(Symmetric(Σ))
    -0.5 * (length(res) * log(2π) + logdet(C) + sum(abs2, C.L \ res))
end

function ais_resid(lm, row)
    run_draw!(lm, row)
    a = sea_level(lm, :ais)
    [a[findfirst(==(y), lm.years)] for y in AY] .- AOBS
end
off(row) = (r = copy(DataFrame(row)); r.antarctic_temp_threshold .= Inf; r[1, :])
function score(lm, row)
    ron, roff = ais_resid(lm, row), ais_resid(lm, off(row))
    sd, rho = Float64(row["sd_ais"]), Float64(row["rho_ais"])
    lon, loff = logl_ais(ron, sd, rho), logl_ais(roff, sd, rho)
    sens = Dict(Symbol("dlogl_$n") => logl_ais(ron, sd, rho, n) - logl_ais(roff, sd, rho, n) for (n, _, _) in NOISE_VARIANTS[2:end])
    (ron = ron, roff = roff, logl_on = lon, logl_off = loff, dlogl = lon - loff, excess_cm = ron[end] - roff[end],
     sens = sens)
end

mkpath(OUTDIR)
post = posterior()
g, o = mean_forcing(FORCING; years = Y0:Y1)
lm = LadrilloModel(g, o; years = Y0:Y1, lws = :central, tap = false, ref = Ladrillo.FIT_REF)
onset = record_crossings(post; last_year = RE, years = Y0:Y1, gmst = g)

## ---- [RESID] and [W=1]: every draw, on and off ------------------------------------------------------------------
println("[RESID] / [W=1] 10,000 draws on the calibration climate ($FORCING, $Y0-$Y1, Antarctic target $(first(AY))-$(last(AY)))")
rf = CSV.read(RESID_FILE, DataFrame; select = vcat(["draw"], ["ais_$y" for y in AY]))
Int.(rf.draw) == collect(1:nrow(post)) || error("residual file draws are not posterior rows 1:$(nrow(post))")
RES = Matrix{Float64}(rf[:, ["ais_$y" for y in AY]])
rows = DataFrame[]
nbad_resid, nbad_w1 = 0, 0
for k in 1:nrow(post)
    ron = ais_resid(lm, post[k, :])
    isequal(ron, RES[k, :]) || (global nbad_resid += 1)
    if onset[k] <= Y1                      # fires by the end of the run: score it
        s = score(lm, post[k, :])
        push!(rows, DataFrame(climate = "calibration ($FORCING mean)", scenario = "", draw = k, file_row = k,
                              config = "", onset = onset[k], sd_ais = post.sd_ais[k], rho_ais = post.rho_ais[k],
                              excess_ais_end_cm = s.excess_cm, logl_on = s.logl_on, logl_off = s.logl_off,
                              dlogl = s.dlogl, weight = exp(s.dlogl); s.sens...))
    else
        isequal(ron, ais_resid(lm, off(post[k, :]))) || (global nbad_w1 += 1)
    end
end
@printf("  %s [RESID] on-residuals == the Table 5 residual file: %d mismatching draws of %d\n",
        nbad_resid == 0 ? "PASS" : "FAIL", nbad_resid, nrow(post))
@printf("  %s [W=1] non-firing draws: on == off residuals exactly (w = 1): %d mismatches of %d\n",
        nbad_w1 == 0 ? "PASS" : "FAIL", nbad_w1, count(>(Y1), onset))

## ---- the joint arm: crossing pairs on their own configs ---------------------------------------------------------
pp = posterior(rows = PROJ_ROWS)
for scen in Ladrillo.SCENARIOS
    mg, mo = mean_forcing(scen; years = 1850:2300)
    cy, cfg, G, O = forcing_cube(scen)
    a = let rng = MersenneTwister(PAIR_SEED), a = Int[]
        while length(a) < nrow(pp); append!(a, randperm(rng, length(cfg))); end; a[1:nrow(pp)]
    end
    iend = findfirst(==(Y1), cy)
    for d in 1:nrow(pp)
        gj = splice(G[:, a[d]], Float64.(mg), cy)[1:iend]
        oj = fastdyn_onset_year(pp[d, :], gj, Y0:Y1)
        oj <= Y1 || continue
        lmj = LadrilloModel(gj, splice(O[:, a[d]], Float64.(mo), cy)[1:iend]; years = Y0:Y1, lws = :central,
                            tap = false, ref = Ladrillo.FIT_REF)
        s = score(lmj, pp[d, :])
        push!(rows, DataFrame(climate = "joint (own spliced config)", scenario = scen, draw = d,
                              file_row = PROJ_ROWS[d], config = cfg[a[d]], onset = oj, sd_ais = pp.sd_ais[d],
                              rho_ais = pp.rho_ais[d], excess_ais_end_cm = s.excess_cm, logl_on = s.logl_on,
                              logl_off = s.logl_off, dlogl = s.dlogl, weight = exp(s.dlogl); s.sens...))
    end
end

df = vcat(rows...)
df.provenance .= "Ladrillo $(pkgversion(Ladrillo)) @ $PKG_COMMIT (posterior L27) | scripts/record_likelihood_weights.jl | " *
    "w = L_ais(fast dynamics on) / L_ais(T_crit = +Inf), the calibrator's AR(1) + correlated-band term " *
    "(obs-corr-len $(OBS_CORR_LEN_YR) yr, Antarctic target $(first(AY))-$(last(AY)), re-ref $(Ladrillo.FIT_REF)) with " *
    "each draw's sd_ais / rho_ais; dlogl_diag* = sensitivity variants (band independent, widened to its 2014-18 " *
    "level $(round(BAND_PRE; digits = 3)) cm after 2018, tripled) | hindcast $Y0-$Y1, land water :central, tap off | calibration climate $FORCING " *
    "mean; joint rows on their own config spliced at $(Ladrillo.SPLICE_YEAR), pairing seed $PAIR_SEED, draws $PROJ_ROWS | " *
    "paleo single assignment, seed $(Ladrillo.PALEO_SEED) | julia $(VERSION)"
f = joinpath(OUTDIR, "record_likelihood_weights.csv")
CSV.write(f, df)
println("\n| climate | scenario | draw | file row | onset | excess AIS $RE (cm) | Δ log L (L27) | weight | Δ log L diag / diag_wide / diag_x3 |")
println("|---|---|---|---|---|---|---|---|---|")
for r in eachrow(df)
    @printf("| %s | %s | %d | %d | %d | %.3f | %.2f | %.3g | %.2f / %.2f / %.2f |\n", r.climate, r.scenario, r.draw,
            r.file_row, r.onset, r.excess_ais_end_cm, r.dlogl, r.weight, r.dlogl_diag, r.dlogl_diag_wide, r.dlogl_diag_x3)
end
println("wrote ", f)
exit(nbad_resid == 0 && nbad_w1 == 0 ? 0 : 1)
