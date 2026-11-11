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
# The BRICK 2.0 comparison runs only if its posterior has been fetched (tools/fetch_brick20_posterior.sh), and so does
# that file's row of the data/MANIFEST.csv check (FETCHED_DATA). Without it a fresh clone reports 2 skipped ("Broken").
#
# Tables 4, A1 and A2 (python/tables/) run only if a Python with NumPy and pandas is found: LADRILLO_PYTHON, else
# ~/climate-env/bin/python, else python3 on the PATH. Their numbers are compared exactly on the PINNED runtime
# (requirements.txt's NumPy and pandas, on macOS arm64) and to RTOL_PORTABLE elsewhere, with the count of values that
# are not bit-identical printed.

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
"""data/MANIFEST.csv rows that are fetched, not shipped: checked when present, skipped when absent (a wrong file FAILS)."""
const FETCHED_DATA = (join(Ladrillo.BRICK20_FILE, "/"),)
const COMPS = (:glaciers, :gis, :ais, :te, :lws, :total)
const SUM_PARTS = (:glaciers, :gis, :ais, :te, :lws)
const SUM_TOL_CM = 1e-6                # scripts/project_joint.jl [SUM]: an identity up to re-referencing float noise
const YEARS = 1850:2300
yi(y) = y - first(YEARS) + 1

agree(a::Real, b::Real) = isequal(a, b) || (!REFERENCE_RUNTIME && abs(a - b) <= RTOL_PORTABLE * max(abs(a), abs(b)))
pairing(n, ncfg) = (rng = MersenneTwister(PAIR_SEED); a = Int[];
                    while length(a) < n; append!(a, randperm(rng, ncfg)); end; a[1:n])
runscript(args...) = success(pipeline(`$JULIA --project=$ROOT $(args)`; stdout = stdout, stderr = stderr))

## Tables A1/A2: the ported research scripts (python/tables/), run in this order, and their research references.
## Table 4 runs on the research hindcast bands in test/reference/ (the BRICK 2.0 band needs the fetched posterior to
## remake here), so it checks the scorecard's arithmetic; the bands themselves are checked by the postpred testset.
const TABLE_DIR = joinpath(ROOT, "python", "tables")
const TABLE_SCRIPTS = ("diag_ais_block_pca.py", "ladrillo_table_a2.py", "ladrillo_prior_posterior_table.py")
const TABLE_GATE = "check_tables_identity.py"
const TABLE4_SCRIPT = "ladrillo_vs_brick20_scorecard.py"
const TABLE4_INPUTS = ("--ladrillo=" * joinpath(REF, "research_postpred_components_timeseries.csv"),
                       "--brick20=" * joinpath(REF, "research_postpred_brick20_components_timeseries.csv"),
                       "--targets=" * joinpath(ROOT, "data", "observations", "calibration_targets.csv"))
const TABLE_A2_REF, TABLE_A2_OUT = "research_ladrillo_table_a2_L27_sub10k.md", "ladrillo_table_a2.md"
const TABLE_A1_REF, TABLE_A1_OUT = "research_ladrillo_prior_posterior_L27.csv", "ladrillo_prior_posterior.csv"
const TABLE_RENAME = ("SLOWP" => "SLOWG", "FAST" => "FASTG")   # the package's glacier-block rename, as whole tokens
const PY_PINNED = ("numpy", "pandas")                           # pinned in requirements.txt
"""A research label in the package's names (gic_a_FAST -> gic_a_FASTG; CUT_FASTDYN, gis_fast untouched)."""
maplabel(s::AbstractString) = foldl((t, p) -> replace(t, Regex("(?<![A-Za-z0-9])$(p.first)(?![A-Za-z0-9])") => p.second),
                                    TABLE_RENAME; init = String(s))
maplabel(x) = x
"""A Python that imports NumPy and pandas, or nothing."""
function table_python()
    for c in (get(ENV, "LADRILLO_PYTHON", ""), joinpath(homedir(), "climate-env", "bin", "python"),
              something(Sys.which("python3"), ""))
        (isempty(c) || !(isfile(c) || Sys.which(c) !== nothing)) && continue
        success(pipeline(`$c -c "import numpy, pandas"`; stdout = devnull, stderr = devnull)) && return c
    end
    nothing
end
drop_provenance_comment(s) = join(filter(l -> !occursin(r"^<!--.*-->$", l), split(s, '\n')), '\n')
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
        if r.file in FETCHED_DATA && !isfile(p)
            @info "$(r.file) not fetched (tools/fetch_brick20_posterior.sh): skipping its MANIFEST check"
            @test_skip false
        else
            @test isfile(p) && bytes2hex(open(sha256, p)) == r.sha256
        end
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

@testset "v0.1 options stay reachable and differ only where documented" begin
    ## paleo: the v0.1 joint-arm assignment changes lambda / T_crit and nothing else
    a, b = posterior(rows = PROJ_ROWS), posterior(rows = PROJ_ROWS, paleo = :per_selection, paleo_block = 500)
    moved = (:antarctic_lambda, :antarctic_temp_threshold)
    @test all(isequal(a[!, c], b[!, c]) for c in names(a) if Symbol(c) ∉ moved)
    @test any(a.antarctic_temp_threshold .!= b.antarctic_temp_threshold)
    ## land water: the v0.1 frame step moves land water by exactly the step before the observed record, Antarctica
    ## through its sea-level feedback, and glaciers / Greenland / thermal expansion not at all
    tg = CSV.read(Ladrillo.datapath("observations", "calibration_targets.csv"), DataFrame)
    obs = Dict(Int(r.year) => Float64(r[Ladrillo.LWS_OBS_COL]) for r in eachrow(tg) if !ismissing(r[Ladrillo.LWS_OBS_COL]))
    fy = minimum(keys(obs)); step = obs[fy]
    mg, mo = mean_forcing("ssp585"; years = YEARS)
    l11, l10 = LadrilloModel(mg, mo; years = YEARS), LadrilloModel(mg, mo; years = YEARS, lws_anchor = :obs_at_first_year)
    for d in (1, 1000, 2000)
        run_draw!(l11, a[d, :]); x = Dict(c => sea_level(l11, c) for c in COMPS)
        run_draw!(l10, a[d, :]); y = Dict(c => sea_level(l10, c) for c in COMPS)
        @test all(isequal(x[c], y[c]) for c in (:glaciers, :gis, :te))
        pre = [yi(t) for t in first(YEARS):(fy - 1)]
        @test maximum(abs.((x[:lws][pre] .- y[:lws][pre]) .- step)) < SUM_TOL_CM
        @test x[:ais] != y[:ais]
    end
end

@testset "record-crossing test (v0.3 opt-in conditioning)" begin
    @test runscript(joinpath(@__DIR__, "conditioning", "check_record_crossing.jl"))
end

@testset "Tables 4, A1, A2 (python/tables) == research" begin
    py = table_python()
    if py === nothing
        @info "No Python with NumPy and pandas found (set LADRILLO_PYTHON, or `pip install -r requirements.txt`): " *
              "skipping the Table 4, A1 and A2 comparison"
        @test_skip false
    else
        vcode = "import " * join(PY_PINNED, ", ") * "; print(" * join(("$m.__version__" for m in PY_PINNED), ", ") * ")"
        got = split(readchomp(`$py -c $vcode`))
        req = Dict(Tuple(split(strip(l), "==")) for l in eachline(joinpath(ROOT, "requirements.txt"))
                   if occursin("==", l) && !startswith(strip(l), "#"))
        pins = [get(req, m, "?") for m in PY_PINNED]
        exact = Sys.isapple() && Sys.ARCH === :aarch64 && got == pins
        rtol = exact ? 0.0 : RTOL_PORTABLE
        println("  python $py ($(join(("$m $v" for (m, v) in zip(PY_PINNED, got)), ", "))): numbers ",
                exact ? "EXACT (pinned runtime)" :
                "to rtol $rtol (not the pinned runtime: requirements.txt $(join(("$m $v" for (m, v) in zip(PY_PINNED, pins)), ", ")) on macOS arm64)")
        out = mktempdir()
        for s in TABLE_SCRIPTS
            @test success(pipeline(`$py $(joinpath(TABLE_DIR, s)) --outdir=$out`; stdout = devnull, stderr = stderr))
        end
        @test success(pipeline(`$py $(joinpath(TABLE_DIR, TABLE4_SCRIPT)) --outdir=$out $(TABLE4_INPUTS)`;
                               stdout = devnull, stderr = stderr))
        ## the identity gate itself, against the research references
        @test success(pipeline(`$py $(joinpath(TABLE_DIR, TABLE_GATE)) --refdir=$REF --pkgdir=$out --rtol=$rtol`;
                               stdout = stdout, stderr = stderr))
        ## and an independent reader (CSV.jl, not pandas) of the same two comparisons
        @test drop_provenance_comment(read(joinpath(REF, TABLE_A2_REF), String)) ==
              drop_provenance_comment(read(joinpath(out, TABLE_A2_OUT), String))
        a = select(CSV.read(joinpath(REF, TABLE_A1_REF), DataFrame), Not(:provenance))
        b = select(CSV.read(joinpath(out, TABLE_A1_OUT), DataFrame), Not(:provenance))
        @test names(a) == names(b) && nrow(a) == nrow(b)
        n = nexact = nbad = 0
        for c in names(a), (x, y) in zip(a[!, c], b[!, c])
            if x isa AbstractFloat && y isa AbstractFloat
                n += 1; nexact += isequal(x, y)
                nbad += !(isequal(x, y) || abs(x - y) <= rtol * max(abs(x), abs(y)))
            else
                nbad += !isequal(maplabel(x), y)
            end
        end
        @printf("  Table A1 (CSV.jl): %d rows x %d cols, %d floats, %d bit-identical, %d outside tolerance or not equal\n",
                nrow(a), ncol(a), n, nexact, nbad)
        @test nbad == 0
    end
end

## The calibration (src/calibration.jl). Two of the gate's checks run here; the chain itself (20,000 iterations, about
## 2 min) is test/calibration/gate_chain_prefix.jl.
@testset "calibration: objective == research L27" begin
    C = Ladrillo.Calibration
    p = C.problem()
    ## (a) the log-posterior at each chain's start row == the first row of that production chain
    for seed in C.SEEDS_CANON
        l = readlines(joinpath(REF, "research_chain_L27_seed$(seed)_row1.csv"))
        hdr, row = split(l[1], ","), split(l[2], ",")
        @test maplabel.(String.(hdr[1:end-2])) == C.param_names(p)
        θ0 = C.start_point(p, seed)
        @test θ0 == parse.(Float64, row[1:end-2])
        lp, lpref = C.logposterior(p, θ0), parse(Float64, row[end-1])
        @printf("  calibration seed %d: logposterior(start) %s, production %s\n", seed, repr(lp), repr(lpref))
        @test agree(lp, lpref)
    end
    ## (b) the parameter table == the research calibrator's own dump of the priors it scores (--dump-priors, Table A1's
    ## input): name, component, mu, sigma and bounds of all 50 parameters
    ref = CSV.read(joinpath(ROOT, "data", "constants", "ladrillo_priors_L27.csv"), DataFrame)
    tab = C.parameter_table(p)
    @test tab.name == String.(ref.name)
    @test tab.comp == String.(ref.comp)
    for c in (:mu, :sigma, :lo, :hi)
        @test all(isequal.(tab[!, c], Float64.(ref[!, c])))
    end
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
