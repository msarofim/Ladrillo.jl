# The observational component targets Ladrillo was calibrated against, and their per-year sigma, on the
# hindcast comparison grid. Shared by the posterior predictive check and the information-criterion residuals so that
# both score against the same numbers.
#
#   ais / gis / steric   Frederikse 2020, extended with GRACE-FO (ice) and NOAA (steric)
#   gsic                 Frederikse 2020 (Marzeion 2015 early segment) spliced to GlaMBIE 2019+, RGI-19 seam removed
#   dang (total)         Frederikse 2020 spliced to NOAA STAR altimetry; NOT fitted (out-of-sample)
#   lws                  the observed land-water budget, ADDED to the modelled total (both models are fitted to
#                        GMSL with land water removed, so the total is only comparable once it is added back)
# Every series is in cm relative to FIT_REF (1995-2005) and starts in FIT_START.

const FIT_START = 1900
const HINDCAST_END = 2026

"""Per-year observational sigma as the likelihood builds it: half the reported 90% band over 1.645, floored at 0.05 cm."""
band_sigma(lo, hi) = max((Float64(hi) - Float64(lo)) / (2 * 1.645), 0.05)

"""
    calibration_targets(; years=FIT_START:HINDCAST_END) -> NamedTuple

`years`; `obs[key]` and `sigma[key]` for key in (:ais, :gsic, :gis, :steric, :dang), NaN where a series has no value;
and `lws`, the observed land-water series. The total's sigma adds the land-water band in quadrature to `dang_sig`.
"""
function calibration_targets(; years=FIT_START:HINDCAST_END)
    FY = collect(Int, years)
    tg = CSV.read(datapath("observations", "calibration_targets.csv"), DataFrame)
    gadj = CSV.read(datapath("observations", "calibration_targets_gsic_adj.csv"), DataFrame)
    tgi(y) = findfirst(==(y), tg.year)
    haveobs(col, y) = (r = tgi(y); r !== nothing && !ismissing(tg[r, col]) && !isnan(Float64(tg[r, col])))
    gd = Dict(Int(gadj[i, :year]) => Float64(gadj[i, :gsic_adj]) for i in 1:nrow(gadj))
    obs_of(key, y) = key === :gsic ? get(gd, y, NaN) : (haveobs(key, y) ? Float64(tg[tgi(y), key]) : NaN)
    function sig(key, y)
        r = tgi(y)
        (r === nothing || isnan(obs_of(key, y))) && return NaN
        key === :dang && return sqrt(Float64(tg[r, :dang_sig])^2 + band_sigma(tg[r, :lws_lo], tg[r, :lws_hi])^2)
        return band_sigma(tg[r, Symbol("$(key)_lo")], tg[r, Symbol("$(key)_hi")])
    end
    keys_ = (:ais, :gsic, :gis, :steric, :dang)
    return (years = FY,
            obs = Dict(k => [obs_of(k, y) for y in FY] for k in keys_),
            sigma = Dict(k => [sig(k, y) for y in FY] for k in keys_),
            lws = [haveobs(:lws, y) ? Float64(tg[tgi(y), :lws]) : NaN for y in FY])
end
