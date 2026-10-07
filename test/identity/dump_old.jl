# Reference side of the v1.0 identity gate: run the RESEARCH kernel (SLR-RFF-BRICK/julia/ladrillo_projection.jl, read
# only) on every case in cases.jl and serialize every component x year x draw.
#
#   julia --project=$SLR/julia_v2 test/identity/dump_old.jl $SLR [case ...]
#
# Draws come from the research repo's own L27 subsample with the OLD column names, read through the old kernel's own
# loader (`ladrillo_posterior`), so the only thing this side shares with the package is cases.jl.

using CSV, DataFrames, Serialization, Printf
const SLR = abspath(ARGS[1])
include(joinpath(SLR, "julia", "ladrillo_projection.jl"))
include(joinpath(@__DIR__, "cases.jl"))
mkpath(OUTDIR)

const OBS = joinpath(SLR, "data", "observations")
mean_old(scen, years) = ([_yearmap(joinpath(OBS, "fair_mean_gmst_$(scen).csv"), "gmst_C")[y] for y in years],
                         [_yearmap(joinpath(OBS, "fair_mean_ohc_$(scen).csv"), "ohc_1e22J")[y] for y in years])

"""Old loader, then the case's rows, then the paleo attachment in the case's blocks (ladrillo_posterior attaches over
whatever it was handed, so it is called once per block exactly as the research drivers did)."""
function draws_old(c)
    hdr = String.(propertynames(CSV.read(LADRILLO_POSTERIOR_CSV, DataFrame; limit=0)))
    full = CSV.read(LADRILLO_POSTERIOR_CSV, DataFrame; select = ladrillo_used_cols(:basins2, hdr))[c.rows, :]
    parts = [ladrillo_native_greenland!(full[lo:min(lo + c.block - 1, nrow(full)), :])
             for lo in 1:c.block:nrow(full)]
    vcat(parts...)
end

function run_case(name, c)
    post = draws_old(c)
    n, ny = nrow(post), length(c.years)
    out = Dict(new => Matrix{Float64}(undef, n, ny) for (_, new) in COMPONENT_MAP)
    function run_into!(idx, g, o)
        bf = ladrillo_setup(ssp = "ssp245", y0 = first(c.years), y1 = last(c.years), gis_variant = :basins2,
                            gmst = g, ohc = o, lws = c.lws)
        c.tap && ladrillo_set_tap!(bf)
        for k in idx
            ladrillo_run_draw!(bf, post[k, :])
            for (old, new) in COMPONENT_MAP
                out[new][k, :] = coalesce.(ladrillo_series(bf, old), NaN)
            end
        end
    end
    mg, mo = mean_old(c.scen, c.years)
    if c.kind === :mean
        run_into!(1:n, mg, mo)
    else
        cg = CSV.read(joinpath(OBS, "fair_cube_gmst_$(c.scen)_raw.csv"), DataFrame)
        co = CSV.read(joinpath(OBS, "fair_cube_ohc_$(c.scen)_raw.csv"), DataFrame)
        cfg = [s for s in String.(propertynames(cg)) if startswith(s, "cfg_")]
        assign = pairing(n, length(cfg))
        for j in unique(assign)
            idx = findall(==(j), assign)
            run_into!(idx, splice_ref(Float64.(cg[!, cfg[j]]), mg, c.years),
                           splice_ref(Float64.(co[!, cfg[j]]), mo, c.years))
        end
    end
    serialize(joinpath(OUTDIR, "old_$(name).jls"), out)
    @printf("[old] %-24s %5d draws x %d years written\n", name, n, ny)
end

for name in (length(ARGS) > 1 ? ARGS[2:end] : sort(collect(keys(CASES))))
    run_case(name, CASES[name])
end
