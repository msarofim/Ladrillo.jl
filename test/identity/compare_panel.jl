# Identity gate 2b (fixed-climate panel): scripts/project_fixed.jl against the research panel.
#   julia --project=. test/identity/compare_panel.jl <research outputs dir> [package outdir]
#   research: ssps_components_2300_L27_tap4p69K_V5p64m_tau800_n2_ws.csv   (julia/project_ssps_components_ladrillo.jl)
#   package : fixed_components_ssp126_ssp245_ssp585.csv                    (scripts/project_fixed.jl)
# Every (year, scenario, component) row: gmst, med, p05, p17, p83, p95 and n_finite with `isequal` (1 ulp fails).
# The research file names the scenario `ssp` and has no provenance column; nothing else may differ.
using CSV, DataFrames, Printf
const RES = abspath(ARGS[1])
const PKG = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(dirname(dirname(@__DIR__)), "outputs", "projections")
a = CSV.read(joinpath(RES, "ssps_components_2300_L27_tap4p69K_V5p64m_tau800_n2_ws.csv"), DataFrame)
b = CSV.read(joinpath(PKG, "fixed_components_ssp126_ssp245_ssp585.csv"), DataFrame)
rename!(a, :ssp => :scenario)
const KEY = [:year, :scenario, :component]
const VALS = [:gmst, :med, :p05, :p17, :p83, :p95, :n_finite]
bad = 0
if nrow(a) != nrow(b)
    @printf("  row counts %d vs %d\n", nrow(a), nrow(b)); global bad += 1
else
    j = innerjoin(select(a, KEY..., VALS...), select(b, KEY..., (v => Symbol(v, :_p) for v in VALS)...), on = KEY)
    nrow(j) == nrow(a) || (@printf("  only %d of %d keys matched\n", nrow(j), nrow(a)); global bad += 1)
    for v in VALS
        nd = count(!, isequal.(j[!, v], j[!, Symbol(v, :_p)]))
        global bad += nd
        @printf("  %-9s %6d rows  %d differ  %s\n", v, nrow(j), nd, nd == 0 ? "IDENTICAL" : "FAIL")
    end
end
@printf("[PANEL-IDENTITY] %s\n", bad == 0 ? "PASS" : "FAIL")
exit(bad == 0 ? 0 : 1)
