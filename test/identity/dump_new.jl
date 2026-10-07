# Package side of the v1.0 identity gate: run Ladrillo on every case in cases.jl and serialize every component x year x
# draw.   julia --project=. test/identity/dump_new.jl [case ...]

using Ladrillo, CSV, DataFrames, Serialization, Printf
include(joinpath(@__DIR__, "cases.jl"))
mkpath(OUTDIR)

function run_case(name, c)
    post = posterior(rows = c.rows, paleo_block = c.block)
    n, ny = nrow(post), length(c.years)
    out = Dict(new => Matrix{Float64}(undef, n, ny) for (_, new) in COMPONENT_MAP)
    function run_into!(idx, g, o)
        lm = LadrilloModel(g, o; years = c.years, lws = c.lws, tap = c.tap)
        for k in idx
            run_draw!(lm, post[k, :])
            for (_, new) in COMPONENT_MAP
                out[new][k, :] = sea_level(lm, new)
            end
        end
    end
    mg, mo = mean_forcing(c.scen; years = c.years)
    if c.kind === :mean
        run_into!(1:n, mg, mo)
    else
        yrs, cfg, G, O = forcing_cube(c.scen)
        yrs == collect(c.years) || error("cube years do not match the case")
        assign = pairing(n, length(cfg))
        for j in unique(assign)
            idx = findall(==(j), assign)
            run_into!(idx, splice(G[:, j], mg, yrs), splice(O[:, j], mo, yrs))
        end
    end
    serialize(joinpath(OUTDIR, "new_$(name).jls"), out)
    @printf("[new] %-24s %5d draws x %d years written\n", name, n, ny)
end

for name in (isempty(ARGS) ? sort(collect(keys(CASES))) : ARGS)
    run_case(name, CASES[name])
end
