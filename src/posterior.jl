# The Ladrillo v1.0 posterior (calibration L27) and the per-draw quantities derived from it.
#
# data/posterior/ladrillo_L27_posterior.csv: 10,000 draws pooled from four 2M-iteration chains (second halves, every
# 400th iteration, in chain order: rows 1-2500 are chain 1, ...). 50 columns. NOT columns, because L27 does not sample
# them: the Antarctic fast-dynamics pair (lambda, T_crit; propagated from the paleo prior as JOINT draws) and gamma
# (held at its paleo median). `posterior()` attaches them.
#
# Two columns are carried in sampling coordinates and mapped back here:
#   Greenland slow channel   (gis_slow_ell, gis_slow_w) -> (gis_alpha_s, gis_beta_s)
#   Antarctic precipitation  ais_precip_u               -> log P0  (in apply_draw!)

const POSTERIOR_FILE = ("posterior", "ladrillo_L27_posterior.csv")
const PALEO_DRAWS_FILE = ("constants", "paleo_fastdyn_draws.csv")
const PALEO_MARGINALS_FILE = ("constants", "paleo_dais_marginals.csv")
"""Seed of the fixed RNG that assigns paleo (lambda, T_crit) rows to posterior draws."""
const PALEO_SEED = 20260920
"""Draws per chain in the 10k posterior file."""
const DRAWS_PER_CHAIN = 2500
const N_CHAINS = 4

"""Greenland slow-channel driver temperature: the 2015-2024 mean of the regional driver (the calibrator's GIS_TBAR,
asserted there at 1.963 K). The slow channel is sampled as r_s = alpha_s*TBAR + beta_s, w = alpha_s*TBAR/r_s."""
const GIS_TBAR_WINDOW = (2015, 2024)
const GIS_ZONE = "south"

function gis_tbar()
    get!(_CACHE, :gis_tbar) do
        tgz = CSV.read(datapath("drivers", "greenland_zone_temperature.csv"), DataFrame)
        t = mean(Float64(tgz[i, GIS_ZONE]) for i in 1:nrow(tgz)
                 if GIS_TBAR_WINDOW[1] <= Int(tgz[i, :year]) <= GIS_TBAR_WINDOW[2])
        abs(t - 1.963) < 5e-3 || error("GIS TBAR = $t disagrees with the calibrator's asserted 1.963 K")
        t
    end
end

"""
    attach_paleo!(df; block=nrow(df))

Attach lambda and T_crit as JOINT paleo draws (rows of paleo_fastdyn_draws.csv chosen by a fixed-seed RNG) and gamma
at its paleo median. The assignment is made in consecutive blocks of `block` rows, the RNG re-seeded for each block.

`block = nrow(df)` (the default) assigns the whole table in one pass. The v1.0 joint-climate arm read its draws chain by
chain and assigned each chain's 500 draws separately (`block = 500`), so its four chains share the same 500 paleo rows;
pass that to reproduce it.
"""
function attach_paleo!(df; block::Int=nrow(df))
    pd = CSV.read(datapath(PALEO_DRAWS_FILE...), DataFrame)
    n = nrow(df)
    lam, tc = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
    for lo in 1:block:n
        hi = min(lo + block - 1, n)
        idx = rand(MersenneTwister(PALEO_SEED), 1:nrow(pd), hi - lo + 1)
        lam[lo:hi] = Float64.(pd.antarctic_lambda[idx])
        tc[lo:hi] = Float64.(pd.antarctic_temp_threshold[idx])
    end
    df.antarctic_lambda = lam
    df.antarctic_temp_threshold = tc
    pm = CSV.read(datapath(PALEO_MARGINALS_FILE...), DataFrame)
    df.antarctic_gamma = fill(Float64(pm.p50[findfirst(==("antarctic_gamma"), pm.param)]), n)
    return df
end

"""Add the native Greenland slow-channel pair from the sampled (gis_slow_ell, gis_slow_w):
r_s = exp(ell), alpha_s = w r_s / TBAR, beta_s = (1 - w) r_s."""
function native_greenland!(df)
    r_s = exp.(Float64.(df.gis_slow_ell)); w_s = Float64.(df.gis_slow_w)
    df.gis_alpha_s = w_s .* r_s ./ gis_tbar()
    df.gis_beta_s  = (1 .- w_s) .* r_s
    return df
end

"""
    posterior(; rows=:, paleo_block=nothing) -> DataFrame

The v1.0 posterior with every per-draw quantity the model needs. `rows` selects draws (e.g. `1:5:10000`, the 2,000
draws of the joint-climate arm) BEFORE the paleo parameters are attached, so the attachment depends only on the
selection. `paleo_block` is passed to `attach_paleo!` (default: one block).
"""
function posterior(; rows=:, paleo_block::Union{Nothing,Int}=nothing)
    df = CSV.read(datapath(POSTERIOR_FILE...), DataFrame)[rows, :]
    attach_paleo!(df; block = paleo_block === nothing ? nrow(df) : paleo_block)
    return native_greenland!(df)
end
