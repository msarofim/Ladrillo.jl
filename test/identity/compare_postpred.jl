# Identity gate 3 (posterior predictive): scripts/posterior_predictive.jl against the research driver's shipped outputs.
#   julia --project=. test/identity/compare_postpred.jl <research outputs dir> [package outdir]
# Compares postpred_L27_{components_timeseries,bias,coverage}.csv (research) with postpred_{...}.csv (package) on every
# column they share except `provenance`, row by row, with `isequal` (so missing == missing and NaN == NaN, and a 1-ulp
# change fails). The package files carry a provenance column the research bias/coverage files lack; that is the only
# column allowed to differ.

using CSV, DataFrames, Printf
const RES = abspath(ARGS[1])
const PKG = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(dirname(dirname(@__DIR__)), "outputs", "hindcast")

bad = 0
for stem in ("components_timeseries", "bias", "coverage")
    a = CSV.read(joinpath(RES, "postpred_L27_$(stem).csv"), DataFrame)
    b = CSV.read(joinpath(PKG, "postpred_$(stem).csv"), DataFrame)
    cols = [c for c in names(a) if c != "provenance"]
    miss = setdiff(cols, names(b))
    extra = setdiff(setdiff(names(b), names(a)), ["provenance"])
    if !isempty(miss) || !isempty(extra) || nrow(a) != nrow(b)
        @printf("  %-22s SHAPE: rows %d vs %d, missing %s, extra %s\n", stem, nrow(a), nrow(b), miss, extra)
        global bad += 1; continue
    end
    nd, w = 0, 0.0
    for c in cols, i in 1:nrow(a)
        x, y = a[i, c], b[i, c]
        isequal(x, y) && continue
        nd += 1
        x isa Real && y isa Real && (w = max(w, abs(Float64(x) - Float64(y))))
    end
    global bad += nd
    @printf("  %-22s %4d rows x %2d cols  %d differ  max|d| %.3e  %s\n", stem, nrow(a), length(cols), nd, w,
            nd == 0 ? "IDENTICAL" : "FAIL")
end
@printf("[POSTPRED-IDENTITY] %s\n", bad == 0 ? "PASS" : "FAIL")
exit(bad == 0 ? 0 : 1)
