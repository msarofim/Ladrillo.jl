# BRICK 2.0 hindcast on the SAME conventions as scripts/posterior_predictive.jl, so a Ladrillo-minus-BRICK reading is
# a module difference and not a driver difference.
#
#   julia --project=. scripts/posterior_predictive_brick20.jl [n_draws] [--outdir=outputs/hindcast]
#
#   run 1850-2026 on the calibration forcing (ssp245harm), saved 1900-2026, re-referenced to FIT_REF (1995-2005)
#   draws: the first n_draws rows of the 10,000-row BRICK 2.0 posterior (default: all)
#   total = ais + glaciers + gis + te + the OBSERVED land-water series. BRICK's own land water is zero before 2019
#           because its calibration target had land water removed, so its components are fitted to GMSL - LWS; adding
#           the observed series back makes the total comparable with the observed total. Not double-counting.
#
# Writes <outdir>/postpred_brick20_components_timeseries.csv: year, {ais,gsic,gis,te,total}_{p5,p50,p95}, provenance.

using Ladrillo, CSV, DataFrames, Printf, Random, Statistics

const Y0, Y1 = 1850, Ladrillo.HINDCAST_END
const FIT_REF = Ladrillo.FIT_REF
const FORCING = "ssp245harm"
const NDRAW_ARG = let p = filter(a -> !startswith(a, "--"), ARGS); isempty(p) ? typemax(Int) : parse(Int, p[1]); end
const OUTDIR = let i = findfirst(a -> startswith(a, "--outdir="), ARGS)
    i === nothing ? joinpath(dirname(@__DIR__), "outputs", "hindcast") : ARGS[i][10:end]
end
mkpath(OUTDIR)
const OUT = joinpath(OUTDIR, "postpred_brick20_components_timeseries.csv")

years = collect(Y0:Y1)
const IB = [findfirst(==(y), years) for y in FIT_REF[1]:FIT_REF[2]]
tgt = Ladrillo.calibration_targets()
FY = tgt.years
myi = [findfirst(==(y), years) for y in FY]
all(isfinite, tgt.lws) || error("observed land water does not cover $(first(FY))-$(last(FY))")

gmst, ohc = mean_forcing(FORCING; years = years)
m = Ladrillo.brick20_model(years = years)
Ladrillo.set_brick20_forcing!(m, gmst, ohc)
post = Ladrillo.brick20_posterior()
ND = min(NDRAW_ARG, nrow(post))
@printf("BRICK 2.0 posterior predictive | %d draws | run %d-%d, saved %d-%d | base %d-%d | forcing %s\n",
        ND, Y0, Y1, first(FY), last(FY), FIT_REF[1], FIT_REF[2], FORCING)

comps = (:ais, :gsic, :gis, :te, :total)
ny = length(FY)
store = Dict(c => Array{Float64}(undef, ND, ny) for c in comps)
t0 = time()
for i in 1:ND
    Ladrillo.apply_brick20!(m, post[i, :])
    run(m)
    ais  = Ladrillo.brick20_sea_level(m, :ais, IB)[myi]
    gsic = Ladrillo.brick20_sea_level(m, :glaciers, IB)[myi]
    gis  = Ladrillo.brick20_sea_level(m, :gis, IB)[myi]
    te   = Ladrillo.brick20_sea_level(m, :te, IB)[myi]
    store[:ais][i, :] = ais; store[:gsic][i, :] = gsic; store[:gis][i, :] = gis; store[:te][i, :] = te
    store[:total][i, :] = ais .+ gsic .+ gis .+ te .+ tgt.lws
end
@printf("  %d draws in %.0fs\n", ND, time() - t0)

band = DataFrame(year = FY)
for c in comps
    band[!, "$(c)_p5"]  = [quantile(store[c][:, j], 0.05) for j in 1:ny]
    band[!, "$(c)_p50"] = [quantile(store[c][:, j], 0.50) for j in 1:ny]
    band[!, "$(c)_p95"] = [quantile(store[c][:, j], 0.95) for j in 1:ny]
end
band[!, "provenance"] = fill("BRICK 2.0 ($(Ladrillo.BRICK20_LABEL)) | Ladrillo $(pkgversion(Ladrillo)) " *
    "scripts/posterior_predictive_brick20.jl | Random.seed!($(Ladrillo.BRICK20_SEED)) immediately before get_model | " *
    "rows 1:$ND | FaIR 2.2.4 (calib 1.6.0) forcing $FORCING | land water = OBSERVED series | run $Y0-$Y1, saved " *
    "$(first(FY))-$(last(FY)), cm rel $(FIT_REF[1])-$(FIT_REF[2]) | julia $(VERSION)", ny)
CSV.write(OUT, band)
println("wrote ", OUT)
