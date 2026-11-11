# Conditioning the propagated paleo draws on the observed record (opt-in; v0.3).
#
# L27 calibrates with lambda / T_crit at their paleo medians, which never reach the fast-dynamics threshold over the
# record. Propagation attaches a JOINT paleo draw to every posterior draw, and a few of those have a threshold low enough
# that fast dynamics fires before the Antarctic record ends, at a rate the record does not show (2.5-8x the same draws
# without it; SLR-RFF-BRICK CHANGELOG 2026-10-08 finding 2). `posterior(; drop_record_crossings=true)` drops them. That
# is rejection sampling, i.e. exact conditioning of the joint (posterior x paleo) draws on "no fast dynamics through
# `record_end()`". It is OFF by default here. The paper's joint-climate arm (v0.3) conditions each draw on its own
# paired climate instead, in scripts/project_joint.jl (see WHOSE CLIMATE below).
#
# The test is the model's own trigger (MimiBRICK antarctic_icesheet: T_ant(t) = (GMST(t-1) - intercept) / coefficient,
# fast dynamics while T_ant > T_crit, never in the first year), evaluated with the arithmetic apply_draw! sets up, so it
# agrees with a model run draw for draw (test/conditioning/check_record_crossing.jl).
#
# WHOSE CLIMATE. The default is the calibration's own hindcast climate, RECORD_FORCING: the climate the likelihood was
# evaluated on, which makes the conditioning a property of the posterior. A run on another climate can cross where that
# one does not. In the joint arm every draw runs on its own spliced FaIR config, and there 3 of the 2,000 projection
# draws cross by 2025 although none does on the mean climate (SLR-RFF-BRICK CHANGELOG 2026-10-08e). To condition a
# joint-climate run, test each draw on the GMST it is actually run with (`fastdyn_onset_year(row, gmst, years)`), AFTER
# pairing: the joint arm pairs draws with configs by position, so dropping draws before pairing re-pairs every later one.

"""The calibration's hindcast climate (the forcing the likelihood was evaluated on): the default climate of the record
test."""
const RECORD_FORCING = "ssp245harm"
"""First year of the record climate: the model's own start year (its first year can never fire)."""
const RECORD_START = 1850

"""Last year of the Antarctic calibration record: the last observed year of the `ais` target in
calibration_targets.csv."""
function record_end()
    get!(_CACHE, :record_end) do
        tg = CSV.read(datapath("observations", "calibration_targets.csv"), DataFrame)
        maximum(Int(tg.year[i]) for i in 1:nrow(tg) if !ismissing(tg.ais[i]) && !isnan(Float64(tg.ais[i])))
    end
end

"""
    fastdyn_onset_year(amp, t_crit, gmst, years) -> Int
    fastdyn_onset_year(row, gmst, years) -> Int

First year of `years` in which Antarctic fast dynamics is on, for a draw with GMST -> T_ant amplification `amp`
(`ais_gmst_amp`) and threshold `t_crit` (`antarctic_temp_threshold`), driven by `gmst` (K rel. 1850-1900, one value per
year of `years`, exactly as fed to `LadrilloModel`). `typemax(Int)` if it never fires within `years`.
"""
function fastdyn_onset_year(amp::Real, t_crit::Real, gmst::AbstractVector{<:Real}, years)
    yrs = collect(Int, years)
    length(gmst) == length(yrs) || error("gmst has $(length(gmst)) values, need $(length(yrs))")
    a, tc = Float64(amp), Float64(t_crit)
    coeff, icept = 1.0 / a, -AIS_TANT0 / a            # as apply_draw! sets ais_temperature_coefficient / _intercept
    for j in 2:length(yrs)
        (Float64(gmst[j-1]) - icept) / coeff > tc && return yrs[j]
    end
    return typemax(Int)
end
fastdyn_onset_year(row, gmst::AbstractVector{<:Real}, years) =
    fastdyn_onset_year(Float64(row["ais_gmst_amp"]), Float64(row["antarctic_temp_threshold"]), gmst, years)

"""The record climate: the mean GMST of `RECORD_FORCING`, RECORD_START through `last_year`."""
record_gmst(last_year::Int=record_end()) = first(mean_forcing(RECORD_FORCING; years = RECORD_START:last_year))

"""
    record_crossings(post; last_year=record_end(), gmst=record_gmst(last_year), years=RECORD_START:last_year)
        -> Vector{Int}

Fast-dynamics onset year of every draw (row) of `post` on one climate, `typemax(Int)` where it does not fire within
`years`. A draw crosses within the record when its onset year is <= `last_year`.
"""
function record_crossings(post::AbstractDataFrame; last_year::Int=record_end(),
                          years=RECORD_START:last_year, gmst::AbstractVector{<:Real}=record_gmst(last(years)))
    last(years) >= last_year || error("years end $(last(years)), before last_year $last_year")
    return [fastdyn_onset_year(post[k, :], gmst, years) for k in 1:nrow(post)]
end
