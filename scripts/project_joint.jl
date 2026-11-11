# Sea-level projections with and without FaIR climate uncertainty, on one scenario.
#
#   julia --project=. scripts/project_joint.jl --scenario=ssp245 [--no-tap] [--outdir=outputs/projections]
#        [--paleo=single|per_selection] [--lws-anchor=zero_at_first_year|obs_at_first_year] [--v0.1]
#
# Two arms on the SAME 2,000 posterior draws:
#   fixed  every draw on the FaIR ensemble-MEAN GMST and ocean heat content
#   joint  each draw paired with ONE of the 841 calibrated FaIR configs (a seeded permutation, cycled), on that
#          config's own path spliced to the mean at 2014 (forcing.jl)
# The joint arm is the paper's reported projection band.
#
# Writes outputs/projections/joint_{cells,draws,paths,gates}_<scenario>[_notap][_paleoper_selection][_lwsobs_at_first_year].csv:
#   cells  median / mean / p05 / p95 per component, horizon and arm
#   draws  every draw's value at each horizon (any statistic can be recomputed without a re-run)
#   paths  median / p05 / p95 per year from 1990
#   gates  the internal checks below
# Every file carries a `provenance` column.
#
# DRAWS: posterior rows 1:5:10000 (= every 2,000th iteration of each chain's second half, 500 per chain).
# v0.2 (the defaults): paleo (lambda, T_crit) rows from ONE assignment over all 10,000 draws; land water starts at 0.
# v0.1 (`--v0.1` = `--paleo=per_selection --lws-anchor=obs_at_first_year`): paleo rows assigned per chain (paleo_block = 500), so the four
# chains reused the same 500 rows; land water carried the 1900 frame step (src/landwater.jl). Non-default settings are
# carried in the output names.

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
"--v1.0" in ARGS && error("--v1.0 was renamed --v0.1 in Ladrillo 1.0.0: it reproduces the first internal build "
                             * "(tag v0.1.0). The defaults are the release settings.")
const V01 = "--v0.1" in ARGS
const PALEO = Symbol(arg("--paleo=", V01 ? "per_selection" : String(Ladrillo.PALEO_ASSIGNMENT)))
PALEO in Ladrillo.PALEO_ASSIGNMENTS || error("--paleo must be one of $(Ladrillo.PALEO_ASSIGNMENTS)")
const ANCHOR = Symbol(arg("--lws-anchor=", V01 ? "obs_at_first_year" : String(Ladrillo.LWS_OBS_ANCHOR)))
ANCHOR in Ladrillo.LWS_OBS_ANCHORS || error("--lws-anchor must be one of $(Ladrillo.LWS_OBS_ANCHORS)")
## RECORD CONDITIONING (v0.3; Marcus 2026-10-08 "own-config, drop", the research driver's rule): the joint arm drops
## every draw whose Antarctic fast dynamics fire by record_end() on the GMST it actually runs with (its own spliced
## config). That is exact rejection = the calibration likelihood's own weighting (CHANGELOG 2026-10-08d). The fixed arm
## (mean climate, on which no projection row crosses) keeps every draw. OFF under --v0.1 (v0.1 did not condition) and
## with --no-record-conditioning (the v0.2 joint arm), which writes "_uncond" names.
const RECORD_CONDITION = !V01 && !("--no-record-conditioning" in ARGS)
const SFX = (TAP ? "" : "_notap") * (PALEO === Ladrillo.PALEO_ASSIGNMENT ? "" : "_paleo$(PALEO)") *
            (ANCHOR === Ladrillo.LWS_OBS_ANCHOR ? "" : "_lws$(ANCHOR)") *
            (RECORD_CONDITION || V01 ? "" : "_uncond")

SCENARIO in Ladrillo.SCENARIOS || error("no shipped forcing for '$SCENARIO'; have $(Ladrillo.SCENARIOS)")
mkpath(OUTDIR)
const PROVENANCE = "Ladrillo $(pkgversion(Ladrillo)) (posterior L27) | scripts/project_joint.jl | scenario $SCENARIO | " *
    "FaIR 2.2.4 (calib 1.6.0) CMIP7-basis forcing, spliced at $(Ladrillo.SPLICE_YEAR) | draws rows $(ROWS) " *
    "(paleo $(PALEO === :per_selection ? "per_selection, blocks of $PALEO_BLOCK" : "single assignment over 10000"), seed $(Ladrillo.PALEO_SEED)) | pairing seed $PAIR_SEED | tap $TAP | land water " *
    ":observed ($ANCHOR) | cm rel $(Ladrillo.REF_PERIOD[1])-$(Ladrillo.REF_PERIOD[2]) | record conditioning " *
    "$(RECORD_CONDITION ? "ON (joint draws firing by $(Ladrillo.record_end()) on their own config dropped)" : "OFF") | " *
    "julia $(VERSION)"

post = PALEO === :per_selection ? posterior(rows = ROWS, paleo = :per_selection, paleo_block = PALEO_BLOCK) : posterior(rows = ROWS)
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
    lm = LadrilloModel(g, o; years = YEARS, tap = TAP, lws_anchor = ANCHOR)
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
const ONSET = [Ladrillo.fastdyn_onset_year(post[k, :], splice(G[:, ASSIGN[k]], mg, yrs), YEARS) for k in 1:NDRAW]
const JKEEP = RECORD_CONDITION ? findall(>(Ladrillo.record_end()), ONSET) : collect(1:NDRAW)
const JSTAT = Dict(c => JOINT[c][JKEEP, :] for c in COMPONENTS)    # every STATISTIC of the joint arm uses this

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
## [RECORD-CONDITIONING] what was dropped, and why (onset year per rejected draw).
let ncross = count(<=(Ladrillo.record_end()), ONSET), drop = setdiff(1:NDRAW, JKEEP)
    push!(gates, ("RECORD-CONDITIONING", "record_end", Float64(Ladrillo.record_end()), "measured"))
    push!(gates, ("RECORD-CONDITIONING", "joint_draws_firing_by_record_end", Float64(ncross), "measured"))
    push!(gates, ("RECORD-CONDITIONING", "joint_draws_dropped", Float64(length(drop)), RECORD_CONDITION ? "ON" : "OFF"))
    for k in drop; push!(gates, ("REJECTED", "draw_$(k)_cfg_$(cfg[ASSIGN[k]])", Float64(ONSET[k]), "onset_year")); end
end
gates.provenance .= PROVENANCE
for r in eachrow(gates); @printf("  [%s] %-32s %.3e  %s\n", r.gate, r.key, r.value, r.verdict); end
any(==("FAIL"), gates.verdict) && error("a gate FAILED; outputs not written")

## ---- outputs -----------------------------------------------------------------------------------------------------
cells = DataFrame(scenario = String[], component = String[], horizon = Int[], arm = String[], n_draws = Int[],
                  med_cm = Float64[], mean_cm = Float64[], p05_cm = Float64[], p95_cm = Float64[])
for c in COMPONENTS, H in HORIZONS, (arm, A) in (("fixed", FIXED[c]), ("joint", JSTAT[c]))
    v = A[:, yi(H)]
    push!(cells, (SCENARIO, String(c), H, arm, length(v), median(v), mean(v), quantile(v, 0.05), quantile(v, 0.95)))
end
draws = DataFrame(draw = Int[], config = String[], component = String[], horizon = Int[], arm = String[],
                  value_cm = Float64[])
for c in COMPONENTS, H in HORIZONS, (arm, A) in (("fixed", FIXED[c]), ("joint", JOINT[c])),
    k in (arm == "joint" ? JKEEP : 1:NDRAW)          # rejected joint draws are not written: conditioning DELETES them
    push!(draws, (k, cfg[ASSIGN[k]], String(c), H, arm, A[k, yi(H)]))
end
paths = DataFrame(year = Int[], component = String[], arm = String[], med_cm = Float64[], p05_cm = Float64[],
                  p95_cm = Float64[])
for c in COMPONENTS, (arm, A) in (("fixed", FIXED[c]), ("joint", JSTAT[c])), (i, y) in enumerate(YEARS)
    y < REPORT0 && continue
    v = A[:, i]
    push!(paths, (y, String(c), arm, median(v), quantile(v, 0.05), quantile(v, 0.95)))
end
for (nm, df) in (("cells", cells), ("draws", draws), ("paths", paths), ("gates", gates))
    df.provenance .= PROVENANCE
    CSV.write(joinpath(OUTDIR, "joint_$(nm)_$(SCENARIO)$(SFX).csv"), df)
end
println("wrote ", joinpath(OUTDIR, "joint_{cells,draws,paths,gates}_$(SCENARIO)$(SFX).csv"))
