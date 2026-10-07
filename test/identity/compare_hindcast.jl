# Identity gate 4 (hindcast products): the BRICK 2.0 posterior predictive and the information-criterion residuals
# against the research drivers' shipped outputs.
#   julia --project=. test/identity/compare_hindcast.jl <research outputs dir> [package outdir]
#
#   research                                         package
#   postpred_oldbrick_components_timeseries.csv      postpred_brick20_components_timeseries.csv
#   ic_hindcast_residuals_ladrillo_L27.csv           ic_hindcast_residuals_ladrillo.csv
#   ic_hindcast_residuals_brick20.csv                ic_hindcast_residuals_brick20.csv
#   ic_hindcast_obs_sigma.csv                        ic_hindcast_obs_sigma.csv
# Every shared column except `provenance`, row by row, with `isequal` (missing == missing, NaN == NaN, 1 ulp fails).

using CSV, DataFrames, Printf
const RES = abspath(ARGS[1])
const PKG = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(dirname(dirname(@__DIR__)), "outputs", "hindcast")
const PAIRS = [("postpred_oldbrick_components_timeseries.csv", "postpred_brick20_components_timeseries.csv"),
               ("ic_hindcast_residuals_ladrillo_L27.csv", "ic_hindcast_residuals_ladrillo.csv"),
               ("ic_hindcast_residuals_brick20.csv", "ic_hindcast_residuals_brick20.csv"),
               ("ic_hindcast_obs_sigma.csv", "ic_hindcast_obs_sigma.csv")]

bad = 0
for (fr, fp) in PAIRS
    a = CSV.read(joinpath(RES, fr), DataFrame)
    b = CSV.read(joinpath(PKG, fp), DataFrame)
    cols = [c for c in names(a) if c != "provenance"]
    miss, extra = setdiff(cols, names(b)), setdiff(setdiff(names(b), names(a)), ["provenance"])
    if !isempty(miss) || !isempty(extra) || nrow(a) != nrow(b)
        @printf("  %-45s SHAPE: rows %d vs %d, missing %s, extra %s\n", fp, nrow(a), nrow(b), miss, extra)
        global bad += 1; continue
    end
    nd, w = 0, 0.0
    for c in cols
        x, y = a[!, c], b[!, c]
        for i in eachindex(x)
            isequal(x[i], y[i]) && continue
            nd += 1
            x[i] isa Real && y[i] isa Real && (w = max(w, abs(Float64(x[i]) - Float64(y[i]))))
        end
    end
    global bad += nd
    @printf("  %-45s %6d rows x %4d cols  %d differ  max|d| %.3e  %s\n", fp, nrow(a), length(cols), nd, w,
            nd == 0 ? "IDENTICAL" : "FAIL")
end
@printf("[HINDCAST-IDENTITY] %s\n", bad == 0 ? "PASS" : "FAIL")
exit(bad == 0 ? 0 : 1)
