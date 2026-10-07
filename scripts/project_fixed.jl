# Sea-level components on the FaIR ensemble-MEAN climate (the fixed-climate panel), per scenario.
#
#   julia --project=. scripts/project_fixed.jl [--scenarios=ssp126,ssp245,ssp585] [--no-tap] [--outdir=...]
#
# 2,000 posterior draws (rows 1:5:10000). v1.0 assigns their paleo (lambda, T_crit) rows in ONE block over the 2,000
# (the joint arm, scripts/project_joint.jl, assigns per chain; see CHANGELOG 2026-10-07).
# Writes outputs/projections/fixed_components_<tag>.csv: year, scenario, component, gmst, med, p05, p17, p83, p95,
# n_finite, provenance -- one row per year from 1990, cm rel. 1995-2014.

using Ladrillo, CSV, DataFrames, Printf, Statistics

arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const SCENS = split(arg("--scenarios=", "ssp126,ssp245,ssp585"), ",")
const TAP = !("--no-tap" in ARGS)
const OUTDIR = arg("--outdir=", joinpath(dirname(@__DIR__), "outputs", "projections"))
const ROWS = 1:5:10000
const YEARS = 1850:2300
const REPORT0 = 1990
const COMPONENTS = [:glaciers, :gis, :ais, :te, :lws, :total]
const LABEL = Dict("ssp126" => "SSP1-2.6", "ssp245" => "SSP2-4.5", "ssp585" => "SSP5-8.5")
const OUT = joinpath(OUTDIR, "fixed_components_$(join(SCENS, "_"))$(TAP ? "" : "_notap").csv")
mkpath(OUTDIR)
const PROVENANCE = "Ladrillo $(pkgversion(Ladrillo)) (posterior L27) | scripts/project_fixed.jl | FaIR 2.2.4 (calib " *
    "1.6.0) CMIP7-basis ensemble-mean forcing | draws rows $(ROWS) (paleo one block, seed $(Ladrillo.PALEO_SEED)) | " *
    "tap $TAP | land water :observed ($(Ladrillo.LWS_OBS_ANCHOR)) | cm rel $(Ladrillo.REF_PERIOD[1])-" *
    "$(Ladrillo.REF_PERIOD[2]) | julia $(VERSION)"

post = posterior(rows = ROWS)
out = DataFrame(year = Int[], scenario = String[], component = String[], gmst = Float64[], med = Float64[],
                p05 = Float64[], p17 = Float64[], p83 = Float64[], p95 = Float64[], n_finite = Int[])
for scen in SCENS
    g, o = mean_forcing(scen; years = YEARS)
    lm = LadrilloModel(g, o; years = YEARS, tap = TAP)
    series = Dict(c => Matrix{Float64}(undef, length(YEARS), nrow(post)) for c in COMPONENTS)
    for j in 1:nrow(post)
        run_draw!(lm, post[j, :])
        for c in COMPONENTS; series[c][:, j] = sea_level(lm, c); end
    end
    for c in COMPONENTS, (i, y) in enumerate(YEARS)
        y >= REPORT0 || continue
        v = filter(isfinite, @view series[c][i, :])
        push!(out, (y, get(LABEL, scen, scen), string(c), lm.gmst[i], median(v), quantile(v, 0.05),
                    quantile(v, 0.17), quantile(v, 0.83), quantile(v, 0.95), length(v)))
    end
    r(c, y) = out[(out.year .== y) .& (out.scenario .== get(LABEL, scen, scen)) .& (out.component .== string(c)), :]
    @printf("%-9s total @2100 %6.1f [%5.1f, %6.1f] cm   @2300 %6.1f cm\n", scen, r(:total, 2100).med[1],
            r(:total, 2100).p17[1], r(:total, 2100).p83[1], r(:total, 2300).med[1])
end
out.provenance .= PROVENANCE
CSV.write(OUT, out)
println("wrote ", OUT)
