# Per-draw hindcast residuals of Ladrillo and BRICK 2.0 against ONE common target set: the input to the
# information-criterion comparison (Table 5; python/ic_ladrillo_vs_brick20.py).
#
#   julia --project=. scripts/ic_hindcast_residuals.jl [ndraw_ladrillo=all] [ndraw_brick20=all] [--outdir=outputs/hindcast]
#
# AIC/BIC need the likelihood at a parameter vector and DIC needs it at every draw, so this writes every draw's residual
# (model - obs, cm) per series and year, with the draw's error-model parameters.
#
# Conventions are those of scripts/posterior_predictive.jl and posterior_predictive_brick20.jl, so these ARE the
# residuals behind their bands (gated below):
#   run 1850-2026 on the calibration forcing (ssp245harm), scored 1900-2026 where a target exists, re-referenced to
#   FIT_REF (1995-2005), targets from calibration_targets(); Ladrillo land water :central and tap off
#   glaciers   Ladrillo: gsic_hind + the per-draw F_unch term; BRICK 2.0: its single glacier reservoir. No delta ramp
#              and no d2 discrepancy on either arm: those are Ladrillo-only likelihood devices.
#   total      ice + steric + the OBSERVED land water on both arms. Out-of-sample for both.
#   draws      evenly thinned, the SAME count on both arms by default (all 10,000): a max-over-draws ln L rises with the
#              number of draws, so unequal counts would favour the larger arm. Posterior row indices are recorded.
#
# Writes <outdir>/ic_hindcast_residuals_{ladrillo,brick20}.csv (draw, error-model columns, <series>_<year>) and
# <outdir>/ic_hindcast_obs_sigma.csv (year, <series>_obs, <series>_sigma), each with a provenance column.

using Ladrillo, CSV, DataFrames, Printf, Random, Statistics

const Y0, Y1 = 1850, Ladrillo.HINDCAST_END
const FIT_REF = Ladrillo.FIT_REF
const FORCING = "ssp245harm"
const POSTPRED_ROWS = 1:5:10000           # scripts/posterior_predictive.jl's default draw set
const POS = filter(a -> !startswith(a, "--"), ARGS)
const NDRAW_L = length(POS) >= 1 ? parse(Int, POS[1]) : typemax(Int)
const NDRAW_B = length(POS) >= 2 ? parse(Int, POS[2]) : typemax(Int)
const OUTDIR = let i = findfirst(a -> startswith(a, "--outdir="), ARGS)
    i === nothing ? joinpath(dirname(@__DIR__), "outputs", "hindcast") : ARGS[i][10:end]
end
mkpath(OUTDIR)
const OUT_L = joinpath(OUTDIR, "ic_hindcast_residuals_ladrillo.csv")
const OUT_B = joinpath(OUTDIR, "ic_hindcast_residuals_brick20.csv")
const OUT_S = joinpath(OUTDIR, "ic_hindcast_obs_sigma.csv")
const PKG = "Ladrillo $(pkgversion(Ladrillo)) scripts/ic_hindcast_residuals.jl"

## series key => target key
const SERIES = [(:ais, :ais), (:gsic, :gsic), (:gis, :gis), (:steric, :steric), (:total, :dang)]
tgt = Ladrillo.calibration_targets()
FY = tgt.years; ny = length(FY)
all(isfinite, tgt.lws) || error("observed land water does not cover $(first(FY))-$(last(FY))")
OBS   = Dict(k => tgt.obs[t] for (k, t) in SERIES)
SIGMA = Dict(k => tgt.sigma[t] for (k, t) in SERIES)

sig = DataFrame(year = FY)
for (k, _) in SERIES
    sig[!, "$(k)_obs"] = OBS[k]; sig[!, "$(k)_sigma"] = SIGMA[k]
end
sig[!, "provenance"] = fill("$PKG | targets calibration_targets.csv + calibration_targets_gsic_adj.csv (gsic) | " *
    "sigma = max((hi-lo)/(2*1.645), 0.05) cm, total adds dang_sig and the land-water band in quadrature | " *
    "cm rel $(FIT_REF[1])-$(FIT_REF[2])", ny)
CSV.write(OUT_S, sig)

function write_residuals(path, draws, noise, resid, prov)
    df = DataFrame(draw = draws)
    for (c, v) in noise; df[!, c] = v; end
    for (k, _) in SERIES, (j, y) in enumerate(FY)
        df[!, "$(k)_$(y)"] = resid[k][:, j]
    end
    df[!, "provenance"] = fill(prov, nrow(df))
    CSV.write(path, df)
end
"""Evenly thinned row indices: every cld(n, k)-th row, at most k of them."""
thin(n, k) = (v = collect(1:cld(n, min(k, n)):n); v[1:min(k, length(v))])

## ---------------------------------------------------------------------------------------------------------------
## Ladrillo
## ---------------------------------------------------------------------------------------------------------------
const NPOST = Ladrillo.DRAWS_PER_CHAIN * Ladrillo.N_CHAINS
ldraws = thin(NPOST, NDRAW_L)
post = posterior(rows = ldraws)
gmst, ohc = mean_forcing(FORCING; years = Y0:Y1)
lm = LadrilloModel(gmst, ohc; years = Y0:Y1, lws = :central, tap = false, ref = FIT_REF)
imy = [findfirst(==(y), lm.years) for y in FY]
@printf("Ladrillo %s (posterior L27) | %d of %d draws | %d-%d scored, base %d-%d | forcing %s\n",
        string(pkgversion(Ladrillo)), nrow(post), NPOST, first(FY), last(FY), FIT_REF[1], FIT_REF[2], FORCING)
resL = Dict(k => Array{Float64}(undef, nrow(post), ny) for (k, _) in SERIES)
t0 = time()
for (i, r) in enumerate(eachrow(post))
    run_draw!(lm, r)
    u = Float64(r["gic_u_unch"])
    ais  = sea_level(lm, :ais)[imy]
    gsic = sea_level(lm, :gsic_hind; funch = u)[imy]
    gis  = sea_level(lm, :gis)[imy]
    te   = sea_level(lm, :te)[imy]
    tot  = sea_level(lm, :glaciers; funch = u)[imy] .+ ais .+ gis .+ te .+ tgt.lws
    resL[:ais][i, :] = ais .- OBS[:ais];  resL[:gsic][i, :] = gsic .- OBS[:gsic]
    resL[:gis][i, :] = gis .- OBS[:gis];  resL[:steric][i, :] = te .- OBS[:steric]
    resL[:total][i, :] = tot .- OBS[:total]
    i % 1000 == 0 && (print("."); flush(stdout))
end
@printf("\n  Ladrillo: %d draws in %.0fs\n", nrow(post), time() - t0)
lnoise = [(c, Float64.(post[!, c])) for c in ("sd_ais", "rho_ais", "sd_gsic", "rho_gsic",
                                              "sd_gis", "rho_gis", "sd_steric", "rho_steric")]
write_residuals(OUT_L, ldraws, lnoise, resL,
    "$PKG | Ladrillo $(pkgversion(Ladrillo)) (posterior L27) rows $(nrow(post)) evenly thinned of $NPOST (paleo :$(Ladrillo.PALEO_ASSIGNMENT) over the file, seed " *
    "$(Ladrillo.PALEO_SEED)) | FaIR 2.2.4 (calib 1.6.0) forcing $FORCING ($(Ladrillo.fair_seed_note(FORCING))) | land water :central | tap off | run " *
    "$Y0-$Y1 scored $(first(FY))-$(last(FY)) | cm rel $(FIT_REF[1])-$(FIT_REF[2]) | gsic = gsic_hind + F_unch vs " *
    "the RGI-19-adjusted target, no delta ramp, no d2 | total = glaciers+ais+gis+te+OBSERVED lws | residual = model - obs")
println("wrote ", OUT_L)

## ---------------------------------------------------------------------------------------------------------------
## BRICK 2.0
## ---------------------------------------------------------------------------------------------------------------
years = collect(Y0:Y1)
ib = [findfirst(==(y), years) for y in FIT_REF[1]:FIT_REF[2]]
myi = [findfirst(==(y), years) for y in FY]
m = Ladrillo.brick20_model(years = years)
Ladrillo.set_brick20_forcing!(m, gmst, ohc)
bpost_full = Ladrillo.brick20_posterior()
bdraws = thin(nrow(bpost_full), NDRAW_B)
bpost = bpost_full[bdraws, :]
@printf("BRICK 2.0 | %d of %d draws | seed %d before get_model\n", nrow(bpost), nrow(bpost_full), Ladrillo.BRICK20_SEED)
resB = Dict(k => Array{Float64}(undef, nrow(bpost), ny) for (k, _) in SERIES)
t0 = time()
for i in 1:nrow(bpost)
    Ladrillo.apply_brick20!(m, bpost[i, :])
    run(m)
    ais  = Ladrillo.brick20_sea_level(m, :ais, ib)[myi]
    gsic = Ladrillo.brick20_sea_level(m, :glaciers, ib)[myi]
    gis  = Ladrillo.brick20_sea_level(m, :gis, ib)[myi]
    te   = Ladrillo.brick20_sea_level(m, :te, ib)[myi]
    tot  = ais .+ gsic .+ gis .+ te .+ tgt.lws
    resB[:ais][i, :] = ais .- OBS[:ais];  resB[:gsic][i, :] = gsic .- OBS[:gsic]
    resB[:gis][i, :] = gis .- OBS[:gis];  resB[:steric][i, :] = te .- OBS[:steric]
    resB[:total][i, :] = tot .- OBS[:total]
    i % 1000 == 0 && (print("."); flush(stdout))
end
@printf("\n  BRICK 2.0: %d draws in %.0fs\n", nrow(bpost), time() - t0)
bnoise = [(c, Float64.(bpost[!, c])) for c in ("sd_glaciers", "rho_glaciers", "sd_greenland", "rho_greenland",
                                               "sd_antarctic", "rho_antarctic", "sd_gmsl", "rho_gmsl")]
write_residuals(OUT_B, bdraws, bnoise, resB,
    "$PKG | BRICK 2.0 ($(Ladrillo.BRICK20_LABEL)) rows $(nrow(bpost)) evenly thinned of " *
    "$(nrow(bpost_full)) | Random.seed!($(Ladrillo.BRICK20_SEED)) immediately before get_model | FaIR 2.2.4 " *
    "(calib 1.6.0) forcing $FORCING ($(Ladrillo.fair_seed_note(FORCING))) | run $Y0-$Y1 scored $(first(FY))-$(last(FY)) | cm rel $(FIT_REF[1])-$(FIT_REF[2]) | " *
    "gsic = gsic_sea_level vs the RGI-19-adjusted target | total = ais+gsic+gis+te+OBSERVED lws | residual = model - obs")
println("wrote ", OUT_B)

## ---------------------------------------------------------------------------------------------------------------
## GATE: the per-draw medians must reproduce the two posterior-predictive p50 series, or these are not the residuals
## behind them. Exact (1e-9) on the postpred's own draw set; a different BRICK 2.0 count is a different subset's
## median, bounded at a quarter of the smallest obs sigma. Run both posterior-predictive scripts first.
## ---------------------------------------------------------------------------------------------------------------
ppL = CSV.read(joinpath(OUTDIR, "postpred_components_timeseries.csv"), DataFrame)
ppB = CSV.read(joinpath(OUTDIR, "postpred_brick20_components_timeseries.csv"), DataFrame)
gidxL = let pos = Dict(d => i for (i, d) in enumerate(ldraws)); [get(pos, d, 0) for d in POSTPRED_ROWS] end
exactL = all(>(0), gidxL)
exactL || (gidxL = collect(1:nrow(post)))
@printf("  [gate] Ladrillo gated on %d rows (%s)\n", length(gidxL),
        exactL ? "exactly the postpred's draw set" : "NOT the postpred's draw set: tolerance gate")
medL = Dict(k => [median(resL[k][gidxL, j]) + OBS[k][j] for j in 1:ny] for (k, _) in SERIES)
medB = Dict(k => [median(resB[k][:, j]) + OBS[k][j] for j in 1:ny] for (k, _) in SERIES)
lmap = Dict(:ais => "ais", :gsic => "glaciers", :gis => "gis", :steric => "te", :total => "total")
bmap = Dict(:ais => "ais", :gsic => "gsic", :gis => "gis", :steric => "te", :total => "total")
ok = true
for (k, _) in SERIES
    fin = findall(j -> isfinite(OBS[k][j]), 1:ny)
    dL = maximum(abs.(medL[k][fin] .- Float64.(ppL[fin, "$(lmap[k])_p50"])))
    dB = maximum(abs.(medB[k][fin] .- Float64.(ppB[fin, "$(bmap[k])_p50"])))
    tolL = exactL ? 1e-9 : 0.25 * minimum(filter(isfinite, SIGMA[k]))
    tolB = nrow(bpost) == nrow(bpost_full) ? 1e-9 : 0.25 * minimum(filter(isfinite, SIGMA[k]))
    passL, passB = dL < tolL, dB < tolB
    @printf("  [gate] %-6s median vs postpred p50: Ladrillo max|d| %.2e vs tol %.0e (%s)  BRICK 2.0 %.2e vs tol %.0e (%s)\n",
            k, dL, tolL, passL ? "PASS" : "FAIL", dB, tolB, passB ? "PASS" : "FAIL")
    global ok &= passL && passB
end
ok || error("residual medians do not reproduce the posterior-predictive p50 series")
println("GATE PASS: per-draw residual medians reproduce the posterior-predictive p50 series")
