# Posterior predictive check: the Ladrillo v1.0 posterior run over the calibration window and compared with the
# observational component targets it was fitted to.
#
#   julia --project=. scripts/posterior_predictive.jl [n_draws] [--outdir=outputs/hindcast]
#
# Conventions are the CALIBRATION's, so a bias means what the fit meant:
#   window     1850-2026 on the calibration forcing (fair_mean_{gmst,ohc}_ssp245harm.csv), land water :central,
#              Greenland tap off (it is a projection-side channel; its onset is above every hindcast warming level)
#   baseline   FIT_REF = 1995-2005 (the calibration re-reference; NOT the 1995-2014 projection baseline)
#   glaciers   model scope gsic_hind = SLOWG + FASTG plus the per-draw uncharted-ice term F_unch, against the glacier
#              target with the RGI-19 seam removed (calibration_targets_gsic_adj.csv)
#   total      modelled ice + steric (with the RGI-19 seam) plus the OBSERVED land-water budget, against the
#              Frederikse 2020 / NOAA STAR spliced total. v1.0 does not fit the total: it is out-of-sample.
#
# Two bands per component:
#   parameter   posterior parameter spread alone
#   predictive  + one draw of the calibrated error model (stationary AR(1) with marginal variance sd^2/(1-rho^2), plus
#               independent per-year observational error), i.e. the residual distribution the likelihood assumes.
#               Coverage of THIS band is the goodness-of-fit statement. NaN for a series with no fitted error model.
#
# Writes <outdir>/postpred_{components_timeseries,bias,coverage}.csv, each with a `provenance` column.
# Draws: the posterior evenly thinned to n_draws (default 2,000 = rows 1:5:10000), paleo rows assigned in one block.

using Ladrillo, CSV, DataFrames, Printf, Random, Statistics

const Y0, Y1 = 1850, 2026
const FIT_REF = Ladrillo.FIT_REF
const FORCING = "ssp245harm"
const FIT_START = 1900                     # every target series starts here
const NOISE_SEED = 2026
const NDRAW_ARG = let p = filter(a -> !startswith(a, "--"), ARGS); isempty(p) ? 2000 : parse(Int, p[1]); end
const OUTDIR = let i = findfirst(a -> startswith(a, "--outdir="), ARGS)
    i === nothing ? joinpath(dirname(@__DIR__), "outputs", "hindcast") : ARGS[i][10:end]
end
mkpath(OUTDIR)
const OUT_BANDS = joinpath(OUTDIR, "postpred_components_timeseries.csv")
const OUT_BIAS  = joinpath(OUTDIR, "postpred_bias.csv")
const OUT_COVER = joinpath(OUTDIR, "postpred_coverage.csv")
const DELTA_END = 1960                     # the glacier delta ramp is zero from here on (v1.0 does not sample it)

## (model component, target column, error-model suffix in the posterior)
const SERIES = [(:ais, :ais, "ais"), (:glaciers, :gsic, "gsic"), (:gis, :gis, "gis"),
                (:te, :steric, "steric"), (:total, :dang, "dang")]

## ---------------------------------------------------------------------------------------------------------------
## targets
## ---------------------------------------------------------------------------------------------------------------
tg   = CSV.read(Ladrillo.datapath("observations", "calibration_targets.csv"), DataFrame)
gadj = CSV.read(Ladrillo.datapath("observations", "calibration_targets_gsic_adj.csv"), DataFrame)
tgi(y) = findfirst(==(y), tg.year)
haveobs(col, y) = (r = tgi(y); r !== nothing && !ismissing(tg[r, col]) && !isnan(Float64(tg[r, col])))

const FY = collect(FIT_START:Y1)
gsic_target = let d = Dict(Int(gadj[i, :year]) => Float64(gadj[i, :gsic_adj]) for i in 1:nrow(gadj))
    Dict(y => get(d, y, NaN) for y in FY)
end
obs_of(key, y) = key === :gsic ? gsic_target[y] : (haveobs(key, y) ? Float64(tg[tgi(y), key]) : NaN)
lws_obs = [haveobs(:lws, y) ? Float64(tg[tgi(y), :lws]) : NaN for y in FY]
delta_ramp = [y < DELTA_END ? (DELTA_END - y) / 10.0 : 0.0 for y in FY]

## Per-year observational sigma, as the likelihood builds it: half the reported band over 1.645, floored at 0.05 cm;
## the total also carries the land-water budget error.
band_sigma(lo, hi) = max((Float64(hi) - Float64(lo)) / (2 * 1.645), 0.05)
function obs_sigma(tcol, y)
    r = tgi(y)
    (r === nothing || isnan(obs_of(tcol, y))) && return NaN
    tcol === :dang && return sqrt(Float64(tg[r, :dang_sig])^2 + band_sigma(tg[r, :lws_lo], tg[r, :lws_hi])^2)
    return band_sigma(tg[r, Symbol("$(tcol)_lo")], tg[r, Symbol("$(tcol)_hi")])
end
const OBS_SIGMA = Dict(tcol => [obs_sigma(tcol, y) for y in FY] for (_, tcol, _) in SERIES)

## ---------------------------------------------------------------------------------------------------------------
## run the posterior
## ---------------------------------------------------------------------------------------------------------------
const NPOST = Ladrillo.DRAWS_PER_CHAIN * Ladrillo.N_CHAINS
const ROWS = let step = cld(NPOST, NDRAW_ARG), r = 1:step:NPOST; r[1:min(NDRAW_ARG, length(r))]; end
post = posterior(rows = ROWS)
## Which series v1.0 was fitted to, read off the posterior: a series is fitted iff it carries sd_* and rho_* columns.
const FITTED = Set(k for (k, _, sfx) in SERIES if ("sd_$sfx" in names(post)) && ("rho_$sfx" in names(post)))
const UNFITTED = [k for (k, _, _) in SERIES if !(k in FITTED)]
isempty(UNFITTED) || println("NOTE: no calibrated error model for $(join(UNFITTED, ", ")): reported OUT-OF-SAMPLE, " *
                             "parameter band only, predictive band NaN")

gmst, ohc = mean_forcing(FORCING; years = Y0:Y1)
lm = LadrilloModel(gmst, ohc; years = Y0:Y1, lws = :central, tap = false, ref = FIT_REF)
imy = [findfirst(==(y), lm.years) for y in FY]
ny = length(FY)
@printf("Ladrillo posterior predictive | v1.0 (L27) | %d draws | %d-%d | base %d-%d | forcing %s\n",
        nrow(post), FIT_START, Y1, FIT_REF[1], FIT_REF[2], FORCING)

model = Dict(k => Array{Float64}(undef, nrow(post), ny) for (k, _, _) in SERIES)
pred  = Dict(k => Array{Float64}(undef, nrow(post), ny) for (k, _, _) in SERIES)
obs_corrected = Array{Float64}(undef, nrow(post), ny)      # glacier target + the per-draw delta ramp (zero in v1.0)
rng = MersenneTwister(NOISE_SEED)
function ar1_plus_obs_error!(rng, out, sd, rho, sigma_y)
    out[1] = sd / sqrt(1 - rho^2) * randn(rng)
    @inbounds for t in 2:length(out)
        out[t] = rho * out[t-1] + sd * randn(rng)
    end
    @inbounds for t in eachindex(out)
        out[t] += isnan(sigma_y[t]) ? 0.0 : sigma_y[t] * randn(rng)
    end
    return out
end
noise = Vector{Float64}(undef, ny)
t0 = time()
for (i, r) in enumerate(eachrow(post))
    run_draw!(lm, r)
    u = Float64(r["gic_u_unch"])
    ais  = sea_level(lm, :ais)[imy]
    gsic = sea_level(lm, :gsic_hind; funch = u)[imy]      # hindcast scope: SLOWG + FASTG + F_unch
    gis  = sea_level(lm, :gis)[imy]
    te   = sea_level(lm, :te)[imy]
    tot  = sea_level(lm, :glaciers; funch = u)[imy] .+ ais .+ gis .+ te .+ lws_obs   # + the RGI-19 seam + obs LWS
    model[:ais][i, :] = ais; model[:glaciers][i, :] = gsic
    model[:gis][i, :] = gis; model[:te][i, :] = te; model[:total][i, :] = tot
    δ = "gic_delta" in names(r) ? Float64(r["gic_delta"]) : 0.0
    obs_corrected[i, :] = [obs_of(:gsic, y) for y in FY] .+ δ .* delta_ramp
    for (key, tcol, sfx) in SERIES
        if key in FITTED
            ar1_plus_obs_error!(rng, noise, Float64(r["sd_$sfx"]), Float64(r["rho_$sfx"]), OBS_SIGMA[tcol])
            pred[key][i, :] = model[key][i, :] .+ noise
        else
            pred[key][i, :] .= NaN
        end
    end
    i % 250 == 0 && (print("."); flush(stdout))
end
@printf("\n  %d draws in %.0fs\n", nrow(post), time() - t0)

## ---------------------------------------------------------------------------------------------------------------
## bands + observations
## ---------------------------------------------------------------------------------------------------------------
const PROVENANCE = "Ladrillo $(pkgversion(Ladrillo)) (posterior L27) | scripts/posterior_predictive.jl | draws rows " *
    "$(ROWS) (paleo single block, seed $(Ladrillo.PALEO_SEED)) | AR(1)+obs noise MersenneTwister(NOISE_SEED=$NOISE_SEED) | " *
    "FaIR 2.2.4 (calib 1.6.0) forcing $FORCING | land water :central | tap off | targets calibration_targets.csv + " *
    "calibration_targets_gsic_adj.csv | cm rel $(FIT_REF[1])-$(FIT_REF[2]) | julia $(VERSION)"
qv(x, p) = (v = filter(isfinite, x); isempty(v) ? NaN : quantile(v, p))
bands = DataFrame(year = FY)
for (key, tcol, _) in SERIES
    bands[!, "$(key)_p05"] = [qv(model[key][:, j], 0.05) for j in 1:ny]
    bands[!, "$(key)_p50"] = [qv(model[key][:, j], 0.50) for j in 1:ny]
    bands[!, "$(key)_p95"] = [qv(model[key][:, j], 0.95) for j in 1:ny]
    bands[!, "$(key)_pred_p05"] = [qv(pred[key][:, j], 0.05) for j in 1:ny]
    bands[!, "$(key)_pred_p95"] = [qv(pred[key][:, j], 0.95) for j in 1:ny]
    o = [obs_of(tcol, y) for y in FY]
    bands[!, "$(key)_obs"] = [isnan(v) ? missing : v for v in o]
end
bands[!, "glaciers_obs_delta_corrected"] = [(v = qv(obs_corrected[:, j], 0.50); isnan(v) ? missing : v) for j in 1:ny]
bands[!, "provenance"] = fill(PROVENANCE, nrow(bands))
CSV.write(OUT_BANDS, bands)

## ---------------------------------------------------------------------------------------------------------------
## bias + coverage. in_sample = FALSE marks a series v1.0 was not fitted to (the total).
## ---------------------------------------------------------------------------------------------------------------
bias = DataFrame(component = String[], year = Int[], obs = Float64[], p50 = Float64[],
                 bias = Float64[], in90 = Bool[], in_sample = Bool[])
function report(key, tcol, y)
    j = findfirst(==(y), FY); j === nothing && return
    o = key === :glaciers ? qv(obs_corrected[:, j], 0.50) : obs_of(tcol, y)
    isnan(o) && return
    p05, p50, p95 = qv(model[key][:, j], 0.05), qv(model[key][:, j], 0.50), qv(model[key][:, j], 0.95)
    inb, insamp = p05 <= o <= p95, key in FITTED
    @printf("  %-9s %d  obs %7.2f  p50 %7.2f  bias %+6.2f  %s%s\n", key, y, o, p50, p50 - o,
            inb ? "in 90%" : "OUT", insamp ? "" : "  [OUT-OF-SAMPLE]")
    push!(bias, (string(key), y, o, p50, p50 - o, inb, insamp))
end
println("\nComponent bias (model p50 - obs, cm; glaciers vs the delta-corrected target)")
for y in (1900, 1950, 2000, 2018)
    println(" @$y:"); for (k, t, _) in SERIES; report(k, t, y); end
end
println(" @series end:")
for (k, t, _) in SERIES
    ys = [y for y in FY if isfinite(k === :glaciers ? obs_of(:gsic, y) : obs_of(t, y))]
    isempty(ys) || report(k, t, maximum(ys))
end
println(" AIS plateau (does the fit track the GRACE-FO pause?):")
for y in (2019, 2020, 2022, 2024, 2025); report(:ais, :ais, y); end
bias[!, "provenance"] = fill(PROVENANCE, nrow(bias))
CSV.write(OUT_BIAS, bias)

println("\n90% coverage over the fit window (share of observed years inside the band)")
cov = DataFrame(component = String[], n_years = Int[], coverage_param = Float64[],
                coverage_pred = Float64[], mean_bias = Float64[], in_sample = Bool[])
for (key, tcol, _) in SERIES
    ins, insp, bs = Bool[], Bool[], Float64[]
    for (j, y) in enumerate(FY)
        o = key === :glaciers ? qv(obs_corrected[:, j], 0.50) : obs_of(tcol, y)
        isnan(o) && continue
        push!(ins, qv(model[key][:, j], 0.05) <= o <= qv(model[key][:, j], 0.95))
        push!(insp, key in FITTED && qv(pred[key][:, j], 0.05) <= o <= qv(pred[key][:, j], 0.95))
        push!(bs, qv(model[key][:, j], 0.50) - o)
    end
    isempty(ins) && continue
    @printf("  %-9s %3d years   parameter %5.1f%%   predictive %5.1f%%   mean bias %+.2f cm%s\n",
            key, length(ins), 100 * mean(ins), key in FITTED ? 100 * mean(insp) : NaN, mean(bs),
            key in FITTED ? "" : "   [OUT-OF-SAMPLE: no fitted error model]")
    push!(cov, (string(key), length(ins), mean(ins), key in FITTED ? mean(insp) : NaN, mean(bs), key in FITTED))
end
cov[!, "provenance"] = fill(PROVENANCE, nrow(cov))
CSV.write(OUT_COVER, cov)
println("\nwrote ", join(basename.([OUT_BANDS, OUT_BIAS, OUT_COVER]), ", "), " -> ", OUTDIR)
