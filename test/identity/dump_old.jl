# Reference side of the identity gate: run the RESEARCH kernel (SLR-RFF-BRICK/julia/ladrillo_projection.jl, read
# only) on every case in cases.jl and serialize every component x year x draw.
#
#   julia --project=$SLR/julia_v2 test/identity/dump_old.jl $SLR [case ...]
#
# The research kernel fixes its land-water anchor and paleo assignment when it is loaded (environment variables
# LADRILLO_LWS_OBS_ANCHOR and LADRILLO_PALEO_ASSIGNMENT, defaulting to v1.1), so one process serves the cases that
# match those settings and SKIPS the rest, saying so. Run it twice to cover both halves:
#   LADRILLO_LWS_OBS_ANCHOR=v1_step LADRILLO_PALEO_ASSIGNMENT=v1 julia ... dump_old.jl $SLR     # the v1.0 cases
#   julia ... dump_old.jl $SLR                                                                    # the v1.1 cases
# A research kernel that predates both settings (v1.0, SLR <= 5d3add0) is treated as (:v1_step, :v1).
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

const OLD_ANCHOR = isdefined(Main, :LWS_OBS_ANCHOR) ? LWS_OBS_ANCHOR : :v1_step
const OLD_PALEO  = isdefined(Main, :LADRILLO_PALEO_ASSIGNMENT) ? LADRILLO_PALEO_ASSIGNMENT : :v1
@printf("research kernel: land-water anchor :%s, paleo assignment :%s\n", OLD_ANCHOR, OLD_PALEO)

"""The case's draws through the research kernel's own loaders.
:v1     the case's rows, then the paleo attachment in the case's blocks (`ladrillo_native_greenland!` attaches over
        whatever it is handed, so it is called once per block exactly as the v1.0 research drivers did).
:single the research loader `ladrillo_posterior` on the whole file (it assigns over all of it), then the case's rows."""
function draws_old(c)
    if c.paleo === :single
        return ladrillo_posterior(path = LADRILLO_POSTERIOR_CSV)[c.rows, :]
    end
    hdr = String.(propertynames(CSV.read(LADRILLO_POSTERIOR_CSV, DataFrame; limit=0)))
    full = CSV.read(LADRILLO_POSTERIOR_CSV, DataFrame; select = ladrillo_used_cols(:basins2, hdr))[c.rows, :]
    parts = [ladrillo_native_greenland!(full[lo:min(lo + c.block - 1, nrow(full)), :])
             for lo in 1:c.block:nrow(full)]
    vcat(parts...)
end

function run_case(name, c)
    if (c.anchor, c.paleo) != (OLD_ANCHOR, OLD_PALEO)
        @printf("[old] %-26s SKIPPED: case is (:%s, :%s), this kernel is (:%s, :%s)\n", name, c.anchor, c.paleo,
                OLD_ANCHOR, OLD_PALEO)
        return
    end
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
