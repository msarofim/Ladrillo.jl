# Ladrillo.jl test suite. Runs from a fresh clone: it needs nothing outside this repository.
#
#   julia --project=. test/runtests.jl                   (about 2 min on an Apple M4; the EXACT check, see below)
#   julia --project=. -e 'using Pkg; Pkg.test()'        (the same tests in Pkg's sandbox, in PORTABLE mode)
#
# WHAT IT TESTS AGAINST. The paper's numbers were produced by the research code (SLR-RFF-BRICK), an independent
# implementation of the same model. test/reference/ holds values those research drivers wrote, copied as text by
# tools/build_test_reference.py (sources, commit and sha256 in test/reference/SOURCES.csv). Re-running the same draws
# here and getting the same numbers is the claim that this package reproduces the paper.
#
# EXACTNESS. On the runtime the references were certified on (REFERENCE_RUNTIME: that platform and Julia version, with
# Julia's default bounds checking) every comparison is exact (`==`), which catches a 1-ulp change to a live parameter
# (mutation-tested, CHANGELOG 2026-10-08e). Anywhere else a value may differ in its last bits, so it must agree to
# RTOL_PORTABLE instead, and the number of values that are not bit-identical is printed. That includes Pkg.test(), which
# forces --check-bounds=yes: the changed codegen re-associates some float operations (measured: <= 1.7e-15 relative,
# the coverage counts unchanged). A threshold draw (fast dynamics, the Greenland tap) can in principle turn such a
# difference into a large one; that FAILS, as it should. Portable mode cannot see a 1-ulp mutation; only the exact run can.
#
# The BRICK 2.0 comparison runs only if its posterior has been fetched (tools/fetch_brick20_posterior.sh).

using Ladrillo, CSV, DataFrames, Printf, Random, SHA, Statistics, Test

const ROOT = dirname(@__DIR__)
const REF = joinpath(@__DIR__, "reference")
const REFERENCE_RUNTIME = Sys.isapple() && Sys.ARCH === :aarch64 && VERSION == v"1.12.6" &&
                          Base.JLOptions().check_bounds == 0          # 0 = Julia's default; Pkg.test() sets 1 (yes)
"""Portable relative tolerance: ~451 annual steps x eps(Float64) ~ 1e-13 of accumulated rounding, with 10x headroom."""
const RTOL_PORTABLE = 1e-12
const JULIA = Base.julia_cmd()
const SCENARIOS = ("ssp245", "ssp585")
const PROJ_ROWS = 1:5:10000            # the Ladrillo joint arm's posterior rows (scripts/project_joint.jl)
const PAIR_SEED = 2026                 # scripts/project_joint.jl, scripts/brick20_joint.jl
const BRICK20_NDRAW = 1000             # scripts/brick20_joint.jl
const COMPS = (:glaciers, :gis, :ais, :te, :lws, :total)
const SUM_PARTS = (:glaciers, :gis, :ais, :te, :lws)
const SUM_TOL_CM = 1e-6                # scripts/project_joint.jl [SUM]: an identity up to re-referencing float noise
const YEARS = 1850:2300
yi(y) = y - first(YEARS) + 1

agree(a::Real, b::Real) = isequal(a, b) || (!REFERENCE_RUNTIME && abs(a - b) <= RTOL_PORTABLE * max(abs(a), abs(b)))
pairing(n, ncfg) = (rng = MersenneTwister(PAIR_SEED); a = Int[];
                    while length(a) < n; append!(a, randperm(rng, ncfg)); end; a[1:n])
runscript(args...) = success(pipeline(`$JULIA --project=$ROOT $(args)`; stdout = stdout, stderr = stderr))
println("Ladrillo $(pkgversion(Ladrillo)) | julia $VERSION | $(Sys.MACHINE) | comparisons ",
        REFERENCE_RUNTIME ? "EXACT (reference runtime)" :
        "PORTABLE, to rtol $RTOL_PORTABLE (not the reference runtime: platform, Julia version or --check-bounds)")

"""Compare a dict (draw, arm, component, horizon) => value with a research draws file; return (n, n_exact, n_bad,
max relative difference)."""
function compare_draws(got, ref)
    n = nexact = nbad = 0; wrel = 0.0
    for r in eachrow(ref)
        v = got[(r.draw, r.arm, Symbol(r.component), r.horizon)]
        n += 1; nexact += isequal(v, r.value_cm); nbad += !agree(v, r.value_cm)
        isequal(v, r.value_cm) || (wrel = max(wrel, abs(v - r.value_cm) / max(abs(v), abs(r.value_cm))))
    end
    (n, nexact, nbad, wrel)
end

@testset "Ladrillo.jl" begin

@testset "data files match data/MANIFEST.csv" begin
    man = CSV.read(joinpath(ROOT, "data", "MANIFEST.csv"), DataFrame)
    for r in eachrow(man)
        p = joinpath(ROOT, "data", r.file)
        @test isfile(p) && bytes2hex(open(sha256, p)) == r.sha256
    end
end

@testset "projections == research draws ($s)" for s in SCENARIOS
    ref = CSV.read(joinpath(REF, "research_joint_draws_$s.csv"), DataFrame)
    draws = sort(unique(ref.draw))
    post = posterior(rows = PROJ_ROWS)
    mg, mo = mean_forcing(s; years = YEARS)
    cy, cfg, G, O = forcing_cube(s)
    assign = pairing(nrow(post), length(cfg))
    @test all(r -> r.arm == "fixed" || r.config == cfg[assign[r.draw]], eachrow(ref))   # the pairing itself
    lmf = LadrilloModel(mg, mo; years = YEARS)
    got = Dict{Tuple{Int,String,Symbol,Int},Float64}()
    sumw = splicew = 0.0
    ylast = Ladrillo.SPLICE_YEAR - Ladrillo.GIS_SHAPE_WINDOW ÷ 2 - 1     # scripts/project_joint.jl [SPLICE-IDENTITY]
    for d in draws
        j = assign[d]
        lmj = LadrilloModel(splice(G[:, j], mg, cy), splice(O[:, j], mo, cy); years = YEARS)
        series = Dict{String,Dict{Symbol,Vector{Float64}}}()
        for (arm, lm) in (("fixed", lmf), ("joint", lmj))
            run_draw!(lm, post[d, :])
            series[arm] = Dict(c => sea_level(lm, c) for c in COMPS)
            for c in COMPS, H in unique(ref.horizon)
                got[(d, arm, c, H)] = series[arm][c][yi(H)]
            end
            sumw = max(sumw, maximum(abs.(sum(series[arm][c] for c in SUM_PARTS) .- series[arm][:total])))
        end
        splicew = max(splicew, maximum(maximum(abs.(series["joint"][c][1:yi(ylast)] .- series["fixed"][c][1:yi(ylast)]))
                                       for c in COMPS))
    end
    n, nexact, nbad, wrel = compare_draws(got, ref)
    @printf("  %s: %d draws x 2 arms, %d values, %d bit-identical, %d outside tolerance, max rel diff %.1e\n", s,
            length(draws), n, nexact, nbad, wrel)
    @test nbad == 0
    @test sumw < SUM_TOL_CM                       # the five components sum to the total
    @test splicew == 0.0                          # the arms are identical until post-splice forcing can reach back
end

@testset "posterior predictive == research" begin
    out = mktempdir()
    @test runscript(joinpath(ROOT, "scripts", "posterior_predictive.jl"), "--outdir=$out")
    for stem in ("components_timeseries", "bias", "coverage")
        a = CSV.read(joinpath(REF, "research_postpred_$stem.csv"), DataFrame)
        b = CSV.read(joinpath(out, "postpred_$stem.csv"), DataFrame)
        cols = [c for c in names(a) if c != "provenance"]
        @test nrow(a) == nrow(b) && issubset(cols, names(b))
        nbad = count(!(isequal(x, y) || (x isa Real && y isa Real && agree(x, y)))
                     for c in cols for (x, y) in zip(a[!, c], b[!, c]))
        @printf("  postpred_%s: %d rows x %d cols, %d outside tolerance\n", stem, nrow(a), length(cols), nbad)
        @test nbad == 0
    end
end

@testset "v1.0 options stay reachable and differ only where documented" begin
    ## paleo: the v1.0 joint-arm assignment changes lambda / T_crit and nothing else
    a, b = posterior(rows = PROJ_ROWS), posterior(rows = PROJ_ROWS, paleo = :v1, paleo_block = 500)
    moved = (:antarctic_lambda, :antarctic_temp_threshold)
    @test all(isequal(a[!, c], b[!, c]) for c in names(a) if Symbol(c) ∉ moved)
    @test any(a.antarctic_temp_threshold .!= b.antarctic_temp_threshold)
    ## land water: the v1.0 frame step moves land water by exactly the step before the observed record, Antarctica
    ## through its sea-level feedback, and glaciers / Greenland / thermal expansion not at all
    tg = CSV.read(Ladrillo.datapath("observations", "calibration_targets.csv"), DataFrame)
    obs = Dict(Int(r.year) => Float64(r[Ladrillo.LWS_OBS_COL]) for r in eachrow(tg) if !ismissing(r[Ladrillo.LWS_OBS_COL]))
    fy = minimum(keys(obs)); step = obs[fy]
    mg, mo = mean_forcing("ssp585"; years = YEARS)
    l11, l10 = LadrilloModel(mg, mo; years = YEARS), LadrilloModel(mg, mo; years = YEARS, lws_anchor = :v1_step)
    for d in (1, 1000, 2000)
        run_draw!(l11, a[d, :]); x = Dict(c => sea_level(l11, c) for c in COMPS)
        run_draw!(l10, a[d, :]); y = Dict(c => sea_level(l10, c) for c in COMPS)
        @test all(isequal(x[c], y[c]) for c in (:glaciers, :gis, :te))
        pre = [yi(t) for t in first(YEARS):(fy - 1)]
        @test maximum(abs.((x[:lws][pre] .- y[:lws][pre]) .- step)) < SUM_TOL_CM
        @test x[:ais] != y[:ais]
    end
end

@testset "record-crossing test (v1.2 opt-in conditioning)" begin
    @test runscript(joinpath(@__DIR__, "conditioning", "check_record_crossing.jl"))
end

@testset "BRICK 2.0 == research draws (ssp245)" begin
    if !isfile(Ladrillo.datapath(Ladrillo.BRICK20_FILE...))
        @info "BRICK 2.0 posterior not fetched (tools/fetch_brick20_posterior.sh): skipping its comparison"
        @test_skip false
    else
        ref = CSV.read(joinpath(REF, "research_brick20_joint_draws_ssp245.csv"), DataFrame)
        P = Ladrillo.brick20_posterior()
        rows = collect(1:max(1, nrow(P) ÷ BRICK20_NDRAW):nrow(P))
        mg, mo = mean_forcing("ssp245"; years = YEARS)
        cy, cfg, G, O = forcing_cube("ssp245")
        assign = pairing(length(rows), length(cfg))
        M = Ladrillo.brick20_model(years = collect(YEARS))
        Ladrillo.set_lws!(M, :observed; anchor = Ladrillo.LWS_OBS_ANCHOR)
        IB = [yi(y) for y in 1995:2014]
        reref(v) = 100 .* (v .- sum(v[IB]) / length(IB))
        var = Dict(:ais => (:antarctic_icesheet, :ais_sea_level), :glaciers => (:glaciers_small_icecaps, :gsic_sea_level),
                   :gis => (:greenland_icesheet, :greenland_sea_level), :te => (:thermal_expansion, :te_sea_level),
                   :lws => (:landwater_storage, :lws_sea_level))
        got = Dict{Tuple{Int,String,Symbol,Int},Float64}()
        for k in sort(unique(ref.draw)), (arm, g, o) in (("fixed", mg, mo),
                ("joint", splice(G[:, assign[k]], mg, cy), splice(O[:, assign[k]], mo, cy)))
            Ladrillo.set_brick20_forcing!(M, g, o)
            Ladrillo.apply_brick20!(M, P[rows[k], :]); run(M)
            s = Dict(c => reref(M[var[c]...]) for c in keys(var))
            s[:total] = s[:ais] .+ s[:glaciers] .+ s[:gis] .+ s[:te] .+ s[:lws]
            for c in COMPS, H in unique(ref.horizon); got[(k, arm, c, H)] = s[c][yi(H)]; end
        end
        n, nexact, nbad, wrel = compare_draws(got, ref)
        @printf("  BRICK 2.0 ssp245: %d values, %d bit-identical, %d outside tolerance, max rel diff %.1e\n", n, nexact,
                nbad, wrel)
        @test nbad == 0
    end
end

end
