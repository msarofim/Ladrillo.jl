# Sea-level projections with and without FaIR climate uncertainty, on one scenario.
#
#   julia --project=. scripts/project_joint.jl --scenario=ssp245 [--no-tap] [--outdir=outputs/projections]
#
# Two arms on the SAME 2,000 posterior draws:
#   fixed  every draw on the FaIR ensemble-MEAN GMST and ocean heat content
#   joint  each draw paired with ONE of the 841 calibrated FaIR configs (a seeded permutation, cycled), on that
#          config's own path spliced to the mean at 2014 (forcing.jl)
# The joint arm is the paper's reported projection band.
#
# Writes outputs/projections/joint_{cells,draws,paths,gates}_<scenario>[_notap].csv:
#   cells  median / mean / p05 / p95 per component, horizon and arm
#   draws  every draw's value at each horizon (any statistic can be recomputed without a re-run)
#   paths  median / p05 / p95 per year from 1990
#   gates  the internal checks below
# Every file carries a `provenance` column.
#
# v1.0 DRAWS: posterior rows 1:5:10000 (= every 2,000th iteration of each chain's second half, 500 per chain). The
# paleo (lambda, T_crit) rows are assigned per chain (paleo_block = 500), as the v1.0 runs did.

using Ladrillo, CSV, DataFrames, Printf, Random, Statistics

arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const SCENARIO = arg("--scenario=", "ssp245")
const TAP = !("--no-tap" in ARGS)
const OUTDIR = arg("--outdir=", joinpath(dirname(@__DIR__), "outputs", "projections"))
const ROWS = 1:5:10000
const PALEO_BLOCK = 500
const PAIR_SEED = 2026
const YEARS = 1850:2300
const HORIZONS = [2100, 2150, 2300]
const COMPONENTS = [:glaciers, :gis, :ais, :te, :lws, :total]
const SUM_PARTS = [:glaciers, :gis, :ais, :te, :lws]
const SUM_TOL_CM = 1e-6                  # an identity in BRICK; float noise only
const REPORT0 = 1990
const SFX = TAP ? "" : "_notap"

SCENARIO in Ladrillo.SCENARIOS || error("no shipped forcing for '$SCENARIO'; have $(Ladrillo.SCENARIOS)")
mkpath(OUTDIR)
const PROVENANCE = "Ladrillo $(pkgversion(Ladrillo)) (posterior L27) | scripts/project_joint.jl | scenario $SCENARIO | " *
    "FaIR 2.2.4 (calib 1.6.0) CMIP7-basis forcing, spliced at $(Ladrillo.SPLICE_YEAR) | draws rows $(ROWS) " *
    "(paleo block $PALEO_BLOCK, seed $(Ladrillo.PALEO_SEED)) | pairing seed $PAIR_SEED | tap $TAP | land water " *
    ":observed ($(Ladrillo.LWS_OBS_ANCHOR)) | cm rel $(Ladrillo.REF_PERIOD[1])-$(Ladrillo.REF_PERIOD[2]) | " *
    "julia $(VERSION)"

post = posterior(rows = ROWS, paleo_block = PALEO_BLOCK)
const NDRAW = nrow(post)
mg, mo = mean_forcing(SCENARIO; years = YEARS)
yrs, cfg, G, O = forcing_cube(SCENARIO)
yrs == collect(YEARS) || error("cube years $(first(yrs))-$(last(yrs)) are not $(YEARS)")
const ASSIGN = let rng = MersenneTwister(PAIR_SEED), a = Int[]
    while length(a) < NDRAW; append!(a, randperm(rng, length(cfg))); end
    a[1:NDRAW]
end
@printf("Ladrillo joint-climate projection | %s | %d draws | %d FaIR configs | tap %s\n",
        SCENARIO, NDRAW, length(cfg), TAP)

alloc() = Dict(c => Matrix{Float64}(undef, NDRAW, length(YEARS)) for c in COMPONENTS)
function run_into!(out, idx, g, o)
    lm = LadrilloModel(g, o; years = YEARS, tap = TAP)
    for k in idx
        run_draw!(lm, post[k, :])
        for c in COMPONENTS; out[c][k, :] = sea_level(lm, c); end
    end
end
const FIXED, JOINT = alloc(), alloc()
run_into!(FIXED, 1:NDRAW, mg, mo)
for (n, j) in enumerate(unique(ASSIGN))
    run_into!(JOINT, findall(==(j), ASSIGN), splice(G[:, j], mg, yrs), splice(O[:, j], mo, yrs))
    n % 100 == 0 && (@printf("  %d configs\n", n); flush(stdout))
end
yi(y) = findfirst(==(y), YEARS)

## ---- gates -------------------------------------------------------------------------------------------------------
gates = DataFrame(gate = String[], key = String[], value = Float64[], verdict = String[])
## [PAIRING] the assignment must use min(configs, draws) configs and must not bias the forcing sample at 2300.
let used = length(unique(ASSIGN)), iref = [yi(y) for y in 1995:2014], i2300 = yi(2300)
    d2300(j) = (s = splice(G[:, j], mg, yrs); s[i2300] - mean(s[iref]))
    drawn, all_ = [d2300(j) for j in ASSIGN], [d2300(j) for j in eachindex(cfg)]
    bias = abs(median(drawn) - median(all_))
    push!(gates, ("PAIRING", "configs_used", used, used == min(length(cfg), NDRAW) ? "PASS" : "FAIL"))
    push!(gates, ("PAIRING", "median_dgmst2300_bias_degC", bias, bias < 0.1 ? "PASS" : "CHECK"))
end
## [SPLICE-IDENTITY] the arms must agree exactly until post-splice forcing can reach back: the Greenland shape law's
## centred running mean (half-width GIS_SHAPE_WINDOW/2) is the first mechanism that does.
let ylast = Ladrillo.SPLICE_YEAR - Ladrillo.GIS_SHAPE_WINDOW ÷ 2 - 1, i1 = yi(ylast)
    w = maximum(maximum(abs.(JOINT[c][:, 1:i1] .- FIXED[c][:, 1:i1])) for c in COMPONENTS)
    push!(gates, ("SPLICE-IDENTITY", "max_cm_through_$(ylast)", w, w < SUM_TOL_CM ? "PASS" : "FAIL"))
end
## [SUM] the five components sum to the total, per draw.
for (arm, R) in (("fixed", FIXED), ("joint", JOINT))
    w = maximum(maximum(abs.(sum(R[c][:, yi(H)] for c in SUM_PARTS) .- R[:total][:, yi(H)])) for H in HORIZONS)
    push!(gates, ("SUM", "$(arm)_max_cm", w, w < SUM_TOL_CM ? "PASS" : "FAIL"))
end
gates.provenance .= PROVENANCE
for r in eachrow(gates); @printf("  [%s] %-32s %.3e  %s\n", r.gate, r.key, r.value, r.verdict); end
any(==("FAIL"), gates.verdict) && error("a gate FAILED; outputs not written")

## ---- outputs -----------------------------------------------------------------------------------------------------
cells = DataFrame(scenario = String[], component = String[], horizon = Int[], arm = String[], n_draws = Int[],
                  med_cm = Float64[], mean_cm = Float64[], p05_cm = Float64[], p95_cm = Float64[])
for c in COMPONENTS, H in HORIZONS, (arm, A) in (("fixed", FIXED[c]), ("joint", JOINT[c]))
    v = A[:, yi(H)]
    push!(cells, (SCENARIO, String(c), H, arm, length(v), median(v), mean(v), quantile(v, 0.05), quantile(v, 0.95)))
end
draws = DataFrame(draw = Int[], config = String[], component = String[], horizon = Int[], arm = String[],
                  value_cm = Float64[])
for c in COMPONENTS, H in HORIZONS, (arm, A) in (("fixed", FIXED[c]), ("joint", JOINT[c])), k in 1:NDRAW
    push!(draws, (k, cfg[ASSIGN[k]], String(c), H, arm, A[k, yi(H)]))
end
paths = DataFrame(year = Int[], component = String[], arm = String[], med_cm = Float64[], p05_cm = Float64[],
                  p95_cm = Float64[])
for c in COMPONENTS, (arm, A) in (("fixed", FIXED[c]), ("joint", JOINT[c])), (i, y) in enumerate(YEARS)
    y < REPORT0 && continue
    v = A[:, i]
    push!(paths, (y, String(c), arm, median(v), quantile(v, 0.05), quantile(v, 0.95)))
end
for (nm, df) in (("cells", cells), ("draws", draws), ("paths", paths), ("gates", gates))
    df.provenance .= PROVENANCE
    CSV.write(joinpath(OUTDIR, "joint_$(nm)_$(SCENARIO)$(SFX).csv"), df)
end
println("wrote ", joinpath(OUTDIR, "joint_{cells,draws,paths,gates}_$(SCENARIO)$(SFX).csv"))
