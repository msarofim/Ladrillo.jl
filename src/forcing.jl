# Climate forcing: FaIR 2.2.4 (calib 1.6.0) GMST and ocean heat content, CMIP7 history, as shipped in data/forcing/.
#
#   fair_mean_{gmst,ohc}_<scenario>.csv       ensemble mean over the 841 calibrated configs (the FIXED-climate arm)
#   fair_cube_{gmst,ohc}_<scenario>_raw.csv   all 841 configs (the JOINT-climate arm)
# GMST in K rel. 1850-1900; OHC in 1e22 J. Regenerate with scripts/fair/ (see docs/forcing.md).
#
# THE SPLICE CONVENTION (joint arm). Through SPLICE_YEAR every config is driven by the ensemble mean (so the hindcast
# is the calibrated one); after it, by the config's own path re-levelled to the mean over REF_PERIOD:
#   x_c(t) = mean(t)                                         t <= SPLICE_YEAR
#   x_c(t) = mean_ref + (raw_c(t) - raw_c,ref)                t >  SPLICE_YEAR

const SPLICE_YEAR = 2014
"""Scenarios with shipped forcing. ssp245harm is the calibration forcing (harmonized SSP2-4.5 emissions)."""
const SCENARIOS = ["ssp126", "ssp245", "ssp585", "vvVL", "vvLN", "vvL", "vvML", "vvM", "vvHL", "vvH"]

"""
    fair_seed_note(scenario) -> String

FaIR's internal variability in the shipped forcing for `scenario` (an SSP, a van Vuuren marker or ssp245harm), for
provenance strings, so that every output depending on a random seed records it. The van Vuuren cubes and means were
built with `stochastic_run=False` (tools/forcing/build_fair_cube_vv_v160.py) and carry no FaIR seed. Every other
shipped forcing (the SSP cubes and means, and the calibration forcing ssp245harm) keeps the calib 1.6.0 posterior's
`stochastic_run=True` with `use_seed=True`: config i is seeded from the posterior's `seed` column (841 distinct values).
The label rests on data/forcing/ being byte-identical to those builds (tools/forcing/check_forcing_identity.py; a vv
cube built with stochastic on fails it).
"""
function fair_seed_note(scenario::AbstractString)
    startswith(scenario, "vv") && return "FaIR stochastic_run=false, no FaIR seed"
    return "FaIR stochastic_run=true, each config seeded from the calib 1.6.0 posterior's seed column"
end

"""Ensemble-mean forcing for `scenario` on `years`: (gmst, ohc)."""
function mean_forcing(scenario::AbstractString; years=1850:2300)
    g = _yearmap(datapath("forcing", "fair_mean_gmst_$(scenario).csv"), "gmst_C")
    o = _yearmap(datapath("forcing", "fair_mean_ohc_$(scenario).csv"), "ohc_1e22J")
    return [g[y] for y in years], [o[y] for y in years]
end

"""Per-config forcing cubes for `scenario`: (years, configs, gmst::Matrix, ohc::Matrix), year x config."""
function forcing_cube(scenario::AbstractString)
    cg = CSV.read(datapath("forcing", "fair_cube_gmst_$(scenario)_raw.csv"), DataFrame)
    co = CSV.read(datapath("forcing", "fair_cube_ohc_$(scenario)_raw.csv"), DataFrame)
    cfg = [c for c in String.(propertynames(cg)) if startswith(c, "cfg_")]
    String.(propertynames(co))[2:end] == cfg || error("gmst and ohc cubes disagree on configs")
    Int.(cg.year) == Int.(co.year) || error("gmst and ohc cubes disagree on years")
    return Int.(cg.year), cfg, Matrix{Float64}(cg[!, cfg]), Matrix{Float64}(co[!, cfg])
end

"""Apply the splice convention to one config's raw path, given the ensemble-mean path on the same `years`."""
function splice(raw::AbstractVector{Float64}, mean_path::AbstractVector{Float64}, years::AbstractVector{Int})
    iref = findall(y -> REF_PERIOD[1] <= y <= REF_PERIOD[2], years)
    mref, cref = mean(mean_path[iref]), mean(raw[iref])
    return [y <= SPLICE_YEAR ? mean_path[i] : mref + (raw[i] - cref) for (i, y) in enumerate(years)]
end
