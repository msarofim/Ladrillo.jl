# Identity gate 2 (products): this package's joint-climate driver against the research driver's own outputs.
#   julia --project=. test/identity/compare_products.jl <research outputs dir> <ladrillo|brick20> <scenario> [...]
# Compares every (draw, config, component, horizon, arm) value of
#   ladrillo  research: scope_slr_fairunc_draws_<scen>_spliced_L27_tap4p69K_V5p64m_tau800.csv
#             package : outputs/projections/joint_draws_<scen>.csv
#   brick20   research: scope_slr_fairunc_draws_<scen>_spliced_oldbrick.csv
#             package : outputs/projections/brick20_draws_<scen>.csv
# with `==` on the parsed Float64 (CSV.jl writes the shortest round-trip representation, so equal text <=> equal bits).
# The research driver read its draws from the four raw 2 GB chains; the package reads the 10k subsample. Passing this
# gate is the proof that the switch (decision 2, 2026-10-07) moved nothing.

using CSV, DataFrames, Printf
const RES = abspath(ARGS[1])
const PKG = joinpath(dirname(dirname(@__DIR__)), "outputs", "projections")
const MODEL = ARGS[2]
MODEL in ("ladrillo", "brick20") || error("model must be ladrillo or brick20")
resfile(s) = MODEL == "ladrillo" ? "scope_slr_fairunc_draws_$(s)_spliced_L27_tap4p69K_V5p64m_tau800.csv" :
                                   "scope_slr_fairunc_draws_$(s)_spliced_oldbrick.csv"
pkgfile(s) = MODEL == "ladrillo" ? "joint_draws_$(s).csv" : "brick20_draws_$(s).csv"

bad = 0
for scen in ARGS[3:end]
    a = CSV.read(joinpath(RES, resfile(scen)), DataFrame)
    b = CSV.read(joinpath(PKG, pkgfile(scen)), DataFrame)
    k = [:draw, :config, :component, :horizon, :arm]
    nrow(a) == nrow(b) || (println("  $scen: row counts $(nrow(a)) vs $(nrow(b))"); global bad += 1; continue)
    j = innerjoin(select(a, k..., :value_cm => :old), select(b, k..., :value_cm => :new), on = k)
    nrow(j) == nrow(a) || (println("  $scen: only $(nrow(j)) of $(nrow(a)) keys matched"); global bad += 1; continue)
    nd = count(!, isequal.(j.old, j.new))
    w = nd == 0 ? 0.0 : maximum(abs.(j.old .- j.new))
    global bad += nd
    @printf("  %-8s %-6s %7d values  %d differ  max|d| %.3e cm  %s\n", MODEL, scen, nrow(j), nd, w, nd == 0 ? "IDENTICAL" : "FAIL")
end
@printf("[PRODUCT-IDENTITY] %s\n", bad == 0 ? "PASS" : "FAIL")
exit(bad == 0 ? 0 : 1)
