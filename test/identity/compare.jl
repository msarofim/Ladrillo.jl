# The identity gate: every component, every year, every draw, `==` (NaN == NaN). Exit 1 on any difference.
#   julia --project=. test/identity/compare.jl [case ...]
using Serialization, Printf
include(joinpath(@__DIR__, "cases.jl"))

function compare(name)
    a = deserialize(joinpath(OUTDIR, "old_$(name).jls"))
    b = deserialize(joinpath(OUTDIR, "new_$(name).jls"))
    bad = 0
    for (_, c) in COMPONENT_MAP
        x, y = a[c], b[c]
        size(x) == size(y) || error("$name/$c: shapes $(size(x)) vs $(size(y))")
        same = isequal.(x, y)                      # NaN-aware, and 0.0 == -0.0
        nd = count(!, same)
        w = nd == 0 ? 0.0 : maximum(abs.(x[.!same] .- y[.!same]))
        bad += nd
        @printf("  %-24s %-10s %8d cells  %8d differ  max|d| %.3e  %s\n", name, c, length(x), nd, w,
                nd == 0 ? "IDENTICAL" : "FAIL")
    end
    return bad
end

total = sum(compare(n) for n in (isempty(ARGS) ? sort(collect(keys(CASES))) : ARGS))
@printf("[IDENTITY] %s\n", total == 0 ? "PASS: bit-identical on every case" : "FAIL: $total cells differ")
exit(total == 0 ? 0 : 1)
