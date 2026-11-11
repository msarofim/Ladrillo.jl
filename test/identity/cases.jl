# Shared case definitions for the identity gate (old research kernel vs this package).
# Included by dump_old.jl (run in SLR-RFF-BRICK's julia_v2 environment) and dump_new.jl (run in this package).
# Neither side may define its own draws, forcing or pairing: they come from here.

using Random, Statistics

## The four v0.1 cases (paleo :per_selection in the blocks the v0.1 drivers used, land water :obs_at_first_year) are the REGRESSION half of
## the gate: the package with its v0.1 options must reproduce the v0.1 research kernel. The `_v11` cases run the v0.2
## defaults (paleo :single = one assignment over all 10,000 draws, then the rows; land water :zero_at_first_year).
## `block` applies to paleo :per_selection only. K3 is hindcast-side (:central land water, which the anchor does not touch), so
## K3 and K3_v11 differ only in the paleo assignment; lambda / T_crit do not reach 1850-2026, so they must be equal.
const CASES = Dict(
    # name   => (forcing, kind, rows, paleo, paleo_block, lws anchor, tap, lws, years)
    "K1_ssp585_mean_tap"     => (scen = "ssp585",     kind = :mean,  rows = 1:10000,  paleo = :per_selection, block = 10000,
                                 anchor = :obs_at_first_year, tap = true,  lws = :observed, years = 1850:2300),
    "K2_ssp245_mean_notap"   => (scen = "ssp245",     kind = :mean,  rows = 1:5:10000, paleo = :per_selection, block = 2000,
                                 anchor = :obs_at_first_year, tap = false, lws = :observed, years = 1850:2300),
    "K3_ssp245harm_central"  => (scen = "ssp245harm", kind = :mean,  rows = 1:10:10000, paleo = :per_selection, block = 1000,
                                 anchor = :obs_at_first_year, tap = false, lws = :central,  years = 1850:2026),
    "K4_vvH_joint_tap"       => (scen = "vvH",        kind = :joint, rows = 1:5:10000, paleo = :per_selection, block = 500,
                                 anchor = :obs_at_first_year, tap = true,  lws = :observed, years = 1850:2300),
    "K1_ssp585_mean_tap_v11"    => (scen = "ssp585",     kind = :mean,  rows = 1:10000,   paleo = :single,
                                    block = nothing, anchor = :zero_at_first_year, tap = true,  lws = :observed,
                                    years = 1850:2300),
    "K2_ssp245_mean_notap_v11"  => (scen = "ssp245",     kind = :mean,  rows = 1:5:10000, paleo = :single,
                                    block = nothing, anchor = :zero_at_first_year, tap = false, lws = :observed,
                                    years = 1850:2300),
    "K3_ssp245harm_central_v11" => (scen = "ssp245harm", kind = :mean,  rows = 1:10:10000, paleo = :single,
                                    block = nothing, anchor = :zero_at_first_year, tap = false, lws = :central,
                                    years = 1850:2026),
    "K4_vvH_joint_tap_v11"      => (scen = "vvH",        kind = :joint, rows = 1:5:10000, paleo = :single,
                                    block = nothing, anchor = :zero_at_first_year, tap = true,  lws = :observed,
                                    years = 1850:2300),
)
"""Old-kernel component name => new component name (the glacier block renames)."""
const COMPONENT_MAP = [(:glaciers, :glaciers), (:gsic_r19, :gsic_r19), (:gsic_slowp, :gsic_slowg),
                       (:gsic_fast, :gsic_fastg), (:gsic_hind, :gsic_hind), (:gis, :gis), (:ais, :ais),
                       (:te, :te), (:lws, :lws), (:total, :total)]
const PAIR_SEED = 2026            # the joint arm's draw -> config permutation (scope_slr_fair_uncertainty.jl)

"""The joint arm's draw -> config assignment: a seeded permutation, cycled."""
function pairing(ndraw::Int, ncfg::Int)
    rng = MersenneTwister(PAIR_SEED)
    a = Int[]
    while length(a) < ndraw; append!(a, randperm(rng, ncfg)); end
    a[1:ndraw]
end

"""The 2014 splice (identical arithmetic to scope_slr_fair_uncertainty.jl `convention`)."""
function splice_ref(raw::Vector{Float64}, mean_path::Vector{Float64}, years)
    iref = findall(y -> 1995 <= y <= 2014, years)
    mref, cref = mean(mean_path[iref]), mean(raw[iref])
    [y <= 2014 ? mean_path[i] : mref + (raw[i] - cref) for (i, y) in enumerate(years)]
end

const OUTDIR = get(ENV, "LADRILLO_IDENTITY_DIR", joinpath(@__DIR__, "out"))
