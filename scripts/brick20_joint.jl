# BRICK 2.0 comparison arm: stock MimiBRICK v2.0.0 with its post-PR#93 posterior, on the SAME climate as Ladrillo.
#
#   julia --project=. scripts/brick20_joint.jl --scenario=ssp245 [--outdir=outputs/projections]
#
# Like-for-like with scripts/project_joint.jl: same FaIR forcing, same 2014 splice, same draw -> config pairing seed,
# same land-water series (:observed), same 1995-2014 reference. Two arms: `fixed` (ensemble-mean climate) and `joint`
# (one FaIR config per draw). 1,000 posterior draws (every 10th row of the 10,000-draw subsample).
#
# MimiBRICK v2.0.0 changed `ais_precipitation₀` to log space; the posterior stores it in linear units, so it is
# log-transformed on the way in (the "precip_log" shim). Total = glaciers + Greenland + Antarctica + TE + land water.
#
# Writes outputs/projections/brick20_{cells,draws,paths}_<scenario>.csv with a provenance column.

using Ladrillo, CSV, DataFrames, Mimi, MimiBRICK, Printf, Random, Statistics

arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const SCENARIO = arg("--scenario=", "ssp245")
const OUTDIR = arg("--outdir=", joinpath(dirname(@__DIR__), "outputs", "projections"))
const NDRAW = 1000
const SEED = 2026                 # global RNG before get_model (whose land-water draw is overwritten anyway)
const PAIR_SEED = 2026
const YEARS = collect(1850:2300)
const HORIZONS = [2100, 2150, 2300]
const COMPS = ["glaciers", "gis", "ais", "te", "lws", "total"]
const REPORT0 = 1990
mkpath(OUTDIR)

const POST = CSV.read(Ladrillo.datapath("brick20", "parameters_subsample_brick.csv"), DataFrame)
const STEP = max(1, nrow(POST) ÷ NDRAW)
const ROWS = collect(1:STEP:nrow(POST))
const PROVENANCE = "BRICK 2.0 (MimiBRICK v2.0.0, post-PR#93 posterior) | scripts/brick20_joint.jl | Ladrillo " *
    "$(pkgversion(Ladrillo)) forcing helpers | scenario $SCENARIO | FaIR 2.2.4 (calib 1.6.0) CMIP7-basis forcing, " *
    "spliced at $(Ladrillo.SPLICE_YEAR) | rows 1:$(STEP):$(nrow(POST)) | pairing seed $PAIR_SEED | land water " *
    ":observed ($(Ladrillo.LWS_OBS_ANCHOR)) | cm rel 1995-2014 | julia $(VERSION)"

mg, mo = mean_forcing(SCENARIO; years = YEARS)
yrs, cfg, G, O = forcing_cube(SCENARIO)
yrs == YEARS || error("cube years are not $(first(YEARS))-$(last(YEARS))")
const ASSIGN = let rng = MersenneTwister(PAIR_SEED), a = Int[]
    while length(a) < length(ROWS); append!(a, randperm(rng, length(cfg))); end
    a[1:length(ROWS)]
end
const IB = [findfirst(==(y), YEARS) for y in 1995:2014]
reref(v) = 100 .* (v .- sum(v[IB]) / length(IB))

"""Apply one BRICK posterior row (MimiBRICK's own parameter names; precip log-transformed for v2.0.0)."""
function update_brick20!(m, r)
    update_param!(m, :antarctic_ocean, :anto_α, r.anto_alpha)
    update_param!(m, :antarctic_ocean, :anto_β, r.anto_beta)
    for (col, sym) in ((:antarctic_s0, :ais_sea_level₀), (:antarctic_bed_height0, :ais_bedheight₀),
                       (:antarctic_slope, :ais_slope), (:antarctic_mu, :ais_μ),
                       (:antarctic_runoff_height0, :ais_runoffline_snowheight₀), (:antarctic_c, :ais_c),
                       (:antarctic_kappa, :ais_κ), (:antarctic_nu, :ais_ν), (:antarctic_flow0, :ais_iceflow₀),
                       (:antarctic_gamma, :ais_γ), (:antarctic_alpha, :ais_α),
                       (:antarctic_temp_threshold, :temperature_threshold), (:antarctic_lambda, :λ))
        update_param!(m, :antarctic_icesheet, sym, r[col])
    end
    update_param!(m, :antarctic_icesheet, :ais_precipitation₀, log(r.antarctic_precip0))
    update_param!(m, :glaciers_small_icecaps, :gsic_β₀, r.glaciers_beta0)
    update_param!(m, :glaciers_small_icecaps, :gsic_v₀, r.glaciers_v0)
    update_param!(m, :glaciers_small_icecaps, :gsic_s₀, r.glaciers_s0)
    update_param!(m, :glaciers_small_icecaps, :gsic_n,  r.glaciers_n)
    update_param!(m, :greenland_icesheet, :greenland_a,  r.greenland_a)
    update_param!(m, :greenland_icesheet, :greenland_b,  r.greenland_b)
    update_param!(m, :greenland_icesheet, :greenland_α,  r.greenland_alpha)
    update_param!(m, :greenland_icesheet, :greenland_β,  r.greenland_beta)
    update_param!(m, :greenland_icesheet, :greenland_v₀, r.greenland_v0)
    update_param!(m, :thermal_expansion, :te_α,  r.thermal_alpha)
    update_param!(m, :thermal_expansion, :te_s₀, r.thermal_s0)
end

alloc() = Dict(c => Matrix{Float64}(undef, length(ROWS), length(YEARS)) for c in COMPS)
Random.seed!(SEED)
const M = MimiBRICK.get_model(ssprcp_scenario = Ladrillo.BUILD_SCENARIO, start_year = first(YEARS), end_year = last(YEARS))
Ladrillo.set_lws!(M, :observed)
function run_into!(store, ks, g, o)
    update_param!(M, :model_global_surface_temperature, g)
    update_param!(M, :thermal_expansion, :ocean_heat_interior, o)
    for k in ks
        update_brick20!(M, POST[ROWS[k], :])
        run(M)
        ais  = reref(M[:antarctic_icesheet,     :ais_sea_level])
        gsic = reref(M[:glaciers_small_icecaps, :gsic_sea_level])
        gis  = reref(M[:greenland_icesheet,     :greenland_sea_level])
        te   = reref(M[:thermal_expansion,      :te_sea_level])
        lws  = reref(M[:landwater_storage,      :lws_sea_level])
        store["ais"][k, :] = ais; store["glaciers"][k, :] = gsic; store["gis"][k, :] = gis
        store["te"][k, :] = te;   store["lws"][k, :] = lws
        store["total"][k, :] = ais .+ gsic .+ gis .+ te .+ lws
    end
end
@printf("BRICK 2.0 joint-climate projection | %s | %d draws | %d FaIR configs\n", SCENARIO, length(ROWS), length(cfg))
const FIXED, JOINT = alloc(), alloc()
run_into!(FIXED, 1:length(ROWS), mg, mo)
let groups = Dict{Int,Vector{Int}}()
    for k in eachindex(ROWS); push!(get!(groups, ASSIGN[k], Int[]), k); end
    for (j, ks) in groups; run_into!(JOINT, ks, splice(G[:, j], mg, yrs), splice(O[:, j], mo, yrs)); end
end
yi(y) = findfirst(==(y), YEARS)

cells = DataFrame(scenario = String[], component = String[], horizon = Int[], arm = String[], n_draws = Int[],
                  med_cm = Float64[], mean_cm = Float64[], p05_cm = Float64[], p95_cm = Float64[])
draws = DataFrame(draw = Int[], config = String[], component = String[], horizon = Int[], arm = String[],
                  value_cm = Float64[])
for c in COMPS, H in HORIZONS, (arm, A) in (("fixed", FIXED[c]), ("joint", JOINT[c]))
    v = A[:, yi(H)]
    push!(cells, (SCENARIO, c, H, arm, length(v), median(v), mean(v), quantile(v, 0.05), quantile(v, 0.95)))
    for k in eachindex(ROWS)
        push!(draws, (k, arm == "joint" ? cfg[ASSIGN[k]] : "MEAN", c, H, arm, A[k, yi(H)]))
    end
end
paths = DataFrame(year = Int[], component = String[], arm = String[], med_cm = Float64[], p05_cm = Float64[],
                  p95_cm = Float64[])
for c in COMPS, (arm, A) in (("fixed", FIXED[c]), ("joint", JOINT[c])), (i, y) in enumerate(YEARS)
    y < REPORT0 && continue
    v = A[:, i]
    push!(paths, (y, c, arm, median(v), quantile(v, 0.05), quantile(v, 0.95)))
end
for (nm, df) in (("cells", cells), ("draws", draws), ("paths", paths))
    df.provenance .= PROVENANCE
    CSV.write(joinpath(OUTDIR, "brick20_$(nm)_$(SCENARIO).csv"), df)
end
println("wrote ", joinpath(OUTDIR, "brick20_{cells,draws,paths}_$(SCENARIO).csv"))
