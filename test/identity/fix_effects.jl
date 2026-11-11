# The v0.1 -> v0.2 change, decomposed on the identity-gate dumps (package side, new_<case>.jls).
#   julia --project=. test/identity/fix_effects.jl
#
# Each v0.2 case is compared with its v0.1 twin. The case pairs isolate the two fixes:
#   K1  rows 1:10000 in ONE paleo block: v1 and :single assign the same rows, so K1 -> K1_v11 is the LAND-WATER fix
#       alone. Gated (exact): glaciers, Greenland, TE IDENTICAL; land water (re-referenced to 1995-2014) rises by
#       exactly obs[first year] before the observed record's first year (v0.1 stepped UP there, so in the 1995-2014
#       frame its earlier years sat that much lower) and moves by float noise from it on; total = land water +
#       Antarctica. Antarctica moves through DAIS's sea-level feedback only; its size is gated against PHYSICS: the
#       median |dAIS| / AIS must be within ANALYTIC_TOL of the analytic flux sensitivity (1 + gamma) / depth times the
#       step (SLR-RFF-BRICK CHANGELOG 2026-10-01f: ~1.1 %/m at a 355 m grounding-line depth). The first version of
#       this gate used the 10-01 "<= 0.012 cm" as its bound and FAILED at 0.049 cm: that number was the transient
#       observed-vs-stylized swap, not this permanent offset (CHANGELOG 2026-10-08, finding 1).
#   K3  hindcast-side (:central land water, which the anchor does not touch), so K3 -> K3_v11 is the PALEO fix alone
#       on 1850-2026. Gated: every component except Antarctica and the total IDENTICAL, and every draw whose Antarctica
#       moves must become identical once T_crit = +Inf in BOTH assignments -- i.e. the only channel is a fast-dynamics
#       threshold crossed before 2027 (4 of the 10,000 draws do; CHANGELOG 2026-10-08, finding 2).
#   K2, K4  both fixes (fixed-climate panel rows, joint arm rows). Reported, not gated.
using Serialization, Printf, CSV, DataFrames, Statistics, Ladrillo
include(joinpath(@__DIR__, "cases.jl"))

const GAMMA_PALEO = 2.83                    # DAIS ice-flow exponent, the paleo median the posterior holds (Table A3)
const GL_DEPTH_M = 355.0                    # grounding-line water depth of the analytic pre-check (CHANGELOG 10-01f)
const ANALYTIC_PER_M = (1 + GAMMA_PALEO) / GL_DEPTH_M   # d ln(flux) / d(sea level), per metre
const ANALYTIC_TOL = 3                      # empirical differing from analytic by > 2-3x is investigated, not accepted
const AIS_FLOOR_CM = 10.0                   # ratios only where Antarctica has contributed something
const FLOAT_CM = 1e-9                       # re-referencing noise on a ~1e2 cm series
const UNMOVED = (:glaciers, :gsic_r19, :gsic_slowg, :gsic_fastg, :gsic_hind, :gis, :te)
## the v0.1 step: the observed series at its first year (cm), read as lws_observed_increments reads it
const FIRST_YEAR, STEP_CM = let tg = CSV.read(Ladrillo.datapath("observations", "calibration_targets.csv"), DataFrame)
    obs = Dict(Int(r.year) => Float64(r[Ladrillo.LWS_OBS_COL]) for r in eachrow(tg) if !ismissing(r[Ladrillo.LWS_OBS_COL]))
    fy = minimum(keys(obs)); (fy, obs[fy])
end

load(name) = deserialize(joinpath(OUTDIR, "new_$(name).jls"))
maxabs(x) = isempty(x) ? 0.0 : maximum(abs, x)
fails = String[]
chk(ok, msg) = (@printf("  %-4s %s\n", ok ? "PASS" : "FAIL", msg); ok || push!(fails, msg))

println("[FIX B ALONE] K1_ssp585_mean_tap -> _v11 (10,000 draws, 1850-2300)")
let a = load("K1_ssp585_mean_tap"), b = load("K1_ssp585_mean_tap_v11"), yrs = collect(CASES["K1_ssp585_mean_tap"].years)
    pre, post = yrs .< FIRST_YEAR, yrs .>= FIRST_YEAR
    for c in UNMOVED
        chk(isequal(a[c], b[c]), "$c identical")
    end
    d = b[:lws] .- a[:lws]
    chk(maxabs(d[:, pre] .- STEP_CM) < FLOAT_CM, @sprintf("lws before %d rises by exactly %.4f cm (max dev %.1e)",
                                                        FIRST_YEAR, STEP_CM, maxabs(d[:, pre] .- STEP_CM)))
    chk(maxabs(d[:, post]) < FLOAT_CM, @sprintf("lws from %d unchanged (max |d| %.1e cm)", FIRST_YEAR, maxabs(d[:, post])))
    dais = b[:ais] .- a[:ais]
    chk(any(!=(0), dais), "ais actually moves (the feedback is live)")
    expect = ANALYTIC_PER_M * STEP_CM / 100
    for H in (2100, 2300)
        i = findfirst(==(H), yrs); ok = abs.(a[:ais][:, i]) .> AIS_FLOOR_CM
        r = median(abs.(dais[ok, i]) ./ abs.(a[:ais][ok, i]))
        chk(expect / ANALYTIC_TOL <= r <= expect * ANALYTIC_TOL,
            @sprintf("ais @%d: median |d|/AIS %.2e vs analytic (1+gamma)/depth x step %.2e (within %dx)",
                     H, r, expect, ANALYTIC_TOL))
    end
    dt = (b[:total] .- a[:total]) .- (d .+ dais)
    chk(maxabs(dt) < FLOAT_CM, @sprintf("total moves by lws + ais only (max residual %.1e cm)", maxabs(dt)))
    for H in (2100, 2300)
        i = findfirst(==(H), yrs)
        @printf("        ais @%d: median d %.2e cm, max |d| %.2e cm\n", H, median(dais[:, i]), maxabs(dais[:, i]))
    end
end

println("[FIX A ON THE HINDCAST] K3_ssp245harm_central -> _v11 (1,000 draws, 1850-2026)")
let n1 = "K3_ssp245harm_central", c = CASES[n1], a = load(n1), b = load("$(n1)_v11")
    for (_, comp) in COMPONENT_MAP
        comp in (:ais, :total) || chk(isequal(a[comp], b[comp]), "$comp identical")
    end
    moved = findall(k -> !isequal(a[:ais][k, :], b[:ais][k, :]), 1:size(a[:ais], 1))
    yrs = collect(c.years)
    fy = [yrs[findfirst(j -> !isequal(a[:ais][k, j], b[:ais][k, j]), eachindex(yrs))] for k in moved]
    @printf("        ais moves in %d of %d draws, first in %s\n", length(moved), size(a[:ais], 1),
            isempty(fy) ? "-" : string(minimum(fy)))
    dt = (b[:total] .- a[:total]) .- (b[:ais] .- a[:ais])
    chk(maxabs(dt) < FLOAT_CM, @sprintf("total moves by ais only (max residual %.1e cm)", maxabs(dt)))
    ## the channel: with the threshold removed in BOTH assignments, the moved draws must agree exactly
    p1 = posterior(rows = c.rows, paleo = :per_selection, paleo_block = c.block); ps = posterior(rows = c.rows)
    p1.antarctic_temp_threshold .= Inf; ps.antarctic_temp_threshold .= Inf
    mg, mo = mean_forcing(c.scen; years = c.years)
    lm = LadrilloModel(mg, mo; years = c.years, lws = c.lws, lws_anchor = c.anchor, tap = c.tap)
    same = all(moved) do k
        run_draw!(lm, p1[k, :]); x = sea_level(lm, :ais)
        run_draw!(lm, ps[k, :]); y = sea_level(lm, :ais)
        isequal(x, y)
    end
    chk(same, "with T_crit = +Inf in both assignments the $(length(moved)) moved draw(s) agree exactly (threshold is the only channel)")
end

println("[BOTH FIXES] reported, not gated")
for (v1, v11) in (("K2_ssp245_mean_notap", "K2_ssp245_mean_notap_v11"), ("K4_vvH_joint_tap", "K4_vvH_joint_tap_v11"))
    a, b, yrs = load(v1), load(v11), collect(CASES[v1].years)
    for H in (2100, 2300), c in (:ais, :total)
        i = findfirst(==(H), yrs)
        @printf("  %-22s %-5s @%d  median %8.3f -> %8.3f cm   p95 %8.3f -> %8.3f cm\n", v1, c, H,
                median(a[c][:, i]), median(b[c][:, i]), quantile(a[c][:, i], 0.95), quantile(b[c][:, i], 0.95))
    end
end
@printf("[FIX-EFFECTS] %s\n", isempty(fails) ? "PASS" : "FAIL: $(length(fails))")
exit(isempty(fails) ? 0 : 1)
