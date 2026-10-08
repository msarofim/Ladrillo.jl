# Checks for the record-crossing test and posterior(; drop_record_crossings=true) (src/record_crossing.jl).
#   julia --project=. test/conditioning/check_record_crossing.jl          (~2 min)
#
# [AGREE-RECORD]  fastdyn_onset_year == the model's own trigger (first year disintegration_rate != 0) for all 10,000
#                 draws on the record climate, RECORD_START:record_end()+1 (one year past the record, so a draw that
#                 first fires the year after it ends is in the comparison too).
# [AGREE-JOINT]   the same on the joint arm (2,000 projection draws, each on its own spliced FaIR config) for
#                 JOINT_SCENARIO: the climate where crossings are most common.
# [BOUNDARY]      POWER. Each NON-crossing draw in BOUNDARY_DRAWS gets T_crit set to its own record-maximum T_ant (as
#                 the model computes it): the trigger is strict, so it must NOT fire. Then T_crit = prevfloat of that:
#                 it must fire by record_end(), in the model AND in the test. A 1-ulp change is enough to flip both.
# [INERT]         drop_record_crossings=false is the v1.1 posterior, and on 1:5:10000 the option drops nothing.
# [DELETION]      the conditioned frame is the unconditioned one with exactly the crossing rows deleted (the paleo
#                 assignment of every kept draw is unchanged), and the dropped set is the documented one.
using Ladrillo, DataFrames, Printf, Random

const JOINT_SCENARIO = "ssp245"
const PAIR_SEED = 2026                       # scripts/project_joint.jl
const PROJ_ROWS = 1:5:10000
const BOUNDARY_DRAWS = (1, 2500, 5000, 7500, 10000)
## the documented crossing set on the record climate (CHANGELOG 2026-10-08c): file row => onset year
const EXPECTED = Dict(9713 => 2021, 4588 => 2025, 4618 => 2025, 4182 => 2026)

fails = String[]
chk(ok, msg) = (@printf("  %-4s %s\n", ok ? "PASS" : "FAIL", msg); ok || push!(fails, msg))
function model_onset(lm, row)
    run_draw!(lm, row)
    dr = lm.m[:antarctic_icesheet, :disintegration_rate]
    j = findfirst(x -> !ismissing(x) && x != 0.0, dr)
    j === nothing ? typemax(Int) : lm.years[j]
end

RE = record_end()
println("record_end() = $RE (last observed year of the Antarctic target); record climate $(Ladrillo.RECORD_FORCING)")
post = posterior()
yrs = Ladrillo.RECORD_START:(RE + 1)
g, o = mean_forcing(Ladrillo.RECORD_FORCING; years = yrs)
lm = LadrilloModel(g, o; years = yrs, lws = :central, tap = false)

println("[AGREE-RECORD] 10,000 draws, $(first(yrs))-$(last(yrs))")
onset = record_crossings(post; last_year = RE, years = yrs, gmst = g)
mo = [model_onset(lm, post[k, :]) for k in 1:nrow(post)]
chk(onset == mo, "test == model trigger on every draw ($(count(onset .!= mo)) mismatches)")
got = Dict(k => onset[k] for k in findall(<=(RE + 1), onset))
chk(got == EXPECTED, "draws firing by $(RE + 1): $(sort(collect(got)))")

println("[AGREE-JOINT] $(JOINT_SCENARIO), $(length(PROJ_ROWS)) projection draws on their own configs")
let pp = posterior(rows = PROJ_ROWS), n = nrow(pp), jy = 1850:2300
    mg, mo_ = mean_forcing(JOINT_SCENARIO; years = jy)
    y, cfg, G, O = forcing_cube(JOINT_SCENARIO)
    a = let rng = MersenneTwister(PAIR_SEED), a = Int[]
        while length(a) < n; append!(a, randperm(rng, length(cfg))); end; a[1:n]
    end
    iend = findfirst(==(RE + 1), y)
    bad, nfire = 0, 0
    for j in unique(a)
        gj = splice(G[:, j], Float64.(mg), y)[1:iend]; oj = splice(O[:, j], Float64.(mo_), y)[1:iend]
        lmj = LadrilloModel(gj, oj; years = 1850:(RE + 1), tap = false)
        for k in findall(==(j), a)
            t = fastdyn_onset_year(pp[k, :], gj, 1850:(RE + 1))
            bad += t != model_onset(lmj, pp[k, :]); nfire += t <= RE
        end
    end
    chk(bad == 0, "test == model trigger on every joint draw ($bad mismatches); $nfire fire by $RE")
end

println("[BOUNDARY] T_crit at / 1 ulp below each draw's record-maximum T_ant")
let ry = Ladrillo.RECORD_START:RE
    rg, ro = mean_forcing(Ladrillo.RECORD_FORCING; years = ry)
    lmr = LadrilloModel(rg, ro; years = ry, lws = :central, tap = false)
    for k in BOUNDARY_DRAWS
        r = copy(post[k:k, :])
        onset[k] > RE || (chk(false, "boundary draw $k already crosses; pick another"); continue)
        amp = Float64(r.ais_gmst_amp[1]); coeff, icept = 1.0 / amp, -Ladrillo.AIS_TANT0 / amp
        tmax = maximum((Float64(rg[j-1]) - icept) / coeff for j in 2:length(ry))
        r.antarctic_temp_threshold .= tmax
        at = (fastdyn_onset_year(r[1, :], rg, ry), model_onset(lmr, r[1, :]))
        r.antarctic_temp_threshold .= prevfloat(tmax)
        below = (fastdyn_onset_year(r[1, :], rg, ry), model_onset(lmr, r[1, :]))
        chk(at == (typemax(Int), typemax(Int)) && below[1] == below[2] <= RE,
            @sprintf("draw %5d: T_crit = max T_ant -> no fire (test, model); 1 ulp below -> fires %d (test) / %d (model)",
                     k, below...))
    end
end

println("[INERT] and [DELETION]")
chk(isequal(posterior(drop_record_crossings = false), post), "drop_record_crossings=false is the v1.1 posterior")
chk(isequal(posterior(rows = PROJ_ROWS, drop_record_crossings = true), posterior(rows = PROJ_ROWS)),
    "on $(PROJ_ROWS) the option drops nothing")
for ly in (RE, RE + 1)
    c = posterior(drop_record_crossings = true, record_last_year = ly)
    want = sort([k for (k, y) in EXPECTED if y <= ly])
    keep = setdiff(1:nrow(post), want)
    chk(isequal(c, post[keep, :]), "record_last_year = $ly: deletes exactly file rows $want, every kept row unchanged")
end

@printf("[RECORD-CROSSING] %s\n", isempty(fails) ? "PASS" : "FAIL: $(length(fails))")
exit(isempty(fails) ? 0 : 1)
