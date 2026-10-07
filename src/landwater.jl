# Land-water storage (LWS).
#
# MimiBRICK's land-water component adds `lws_random_sample[t]` (m) per year after `first_projection_year`, starting
# from `lws₀`. Its stock behaviour draws that series fresh and UNSEEDED in get_model; Ladrillo always overwrites it.
#
#   :observed  PROJECTIONS. The Frederikse 2020 + GRACE/GRACE-FO land-water series (the `lws` column of the calibration
#              target) through LWS_OBS_LAST_REAL_YEAR, then LWS_MEAN per year.
#   :central   HINDCAST-SIDE runs (calibration, posterior predictive): a constant LWS_MEAN per year, the convention the
#              likelihood was fitted under.
#
# KNOWN ISSUE IN v1.0 (the land-water frame step). The observed series is in the calibration target's 1995-2005 frame,
# while every other component starts at 0 in 1850, so under :observed the component steps from 0 to obs[first year]
# (+1.565 cm) in 1900 and the Antarctic sea-level feedback sees that step. Measured on all 2,000 L27 projection draws:
# <= 0.012 cm of Antarctic contribution at any year to 2300; every reported total is re-referenced to 1995-2014, so
# nothing else moves. v1.0 ships WITH the step (LWS_OBS_ANCHOR = :v1_step) so its published numbers reproduce; the fix
# (:zero_at_first_year, the component starts at 0 and carries the observed increments) is scheduled for v1.1.

const LWS_MEAN = 0.0003                 # m/yr, MimiBRICK's mean land-water rate
const LWS_OBS_COL = :lws                # cm, rel. 1995-2005 (the calibration target's own frame)
const LWS_OBS_LAST_REAL_YEAR = 2023     # GRACE-FO record end; the target holds 2023 flat to 2026
const LWS_OBS_ANCHOR = :v1_step
const LWS_OBS_ANCHORS = (:v1_step, :zero_at_first_year)
const LWS_MODES = (:observed, :central)

"""
    lws_observed_increments(years; anchor=LWS_OBS_ANCHOR) -> (first_year, lws0_m, increments_m)

The observed land-water series as MimiBRICK's land-water component consumes it: 0 before `first_year`, `lws0_m` at
`first_year`, plus `increments_m[t]` each year after. Refuses a target whose record extends past
LWS_OBS_LAST_REAL_YEAR with real (non-held) values, so a target update cannot silently change the splice.
"""
function lws_observed_increments(years::AbstractVector{Int}; anchor::Symbol=LWS_OBS_ANCHOR)
    anchor in LWS_OBS_ANCHORS || error("anchor must be one of $(LWS_OBS_ANCHORS), got :$anchor")
    tg = CSV.read(datapath("observations", "calibration_targets.csv"), DataFrame)
    obs = Dict(Int(r.year) => Float64(r[LWS_OBS_COL]) for r in eachrow(tg) if !ismissing(r[LWS_OBS_COL]))
    oyrs = sort(collect(keys(obs)))
    first_year = oyrs[1]
    all(y -> haskey(obs, y) && isfinite(obs[y]), first_year:LWS_OBS_LAST_REAL_YEAR) ||
        error("observed land water is not finite over $(first_year)-$(LWS_OBS_LAST_REAL_YEAR)")
    held = [y for y in (LWS_OBS_LAST_REAL_YEAR+1):oyrs[end] if obs[y] != obs[LWS_OBS_LAST_REAL_YEAR]]
    isempty(held) || error("the target carries non-held values after $(LWS_OBS_LAST_REAL_YEAR) ($(held)); " *
                           "the record was extended -- move LWS_OBS_LAST_REAL_YEAR deliberately")
    first_year in years || error("model years $(first(years))-$(last(years)) do not contain $(first_year)")
    v = zeros(length(years))
    for (i, y) in enumerate(years)
        if y > first_year && y <= LWS_OBS_LAST_REAL_YEAR
            v[i] = (obs[y] - obs[y-1]) / 100.0
        elseif y > LWS_OBS_LAST_REAL_YEAR
            v[i] = LWS_MEAN
        end
    end
    return first_year, (anchor === :v1_step ? obs[first_year] / 100.0 : 0.0), v
end

"""Set the land-water series of a built model `m` to mode `lws` (:observed or :central)."""
function set_lws!(m, lws::Symbol; anchor::Symbol=LWS_OBS_ANCHOR)
    lws in LWS_MODES || error("lws must be one of $(LWS_MODES), got :$lws")
    years = Int.(Mimi.dim_keys(m, :time))
    if lws === :observed
        first_year, lws0_m, v = lws_observed_increments(years; anchor=anchor)
        update_param!(m, :landwater_storage, :first_projection_year, first_year)
        update_param!(m, :landwater_storage, :lws₀, lws0_m)
        update_param!(m, :landwater_storage, :lws_random_sample, v)
    else
        update_param!(m, :landwater_storage, :lws_random_sample, fill(LWS_MEAN, length(years)))
    end
    return m
end
