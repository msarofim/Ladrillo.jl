# The BRICK 2.0 comparison arm: stock MimiBRICK v2.0.0 with the BRICK parameters of its published joint SNEASY-BRICK
# posterior subsample (Zenodo 10.5281/zenodo.20592337, `parameters_subsample_sneasybrick.csv`, 10,000 draws; Marcus
# 2026-10-15). The 16 SNEASY columns are present in the file and not read: the arm is driven by FaIR.
#
# The posterior is NOT redistributed with this package: tools/fetch_brick20_posterior.sh fetches it into
# data/brick20/ and checks its sha256.
#
# The file stores `antarctic_precip0` in log space, which is what MimiBRICK v2.0.0's `ais_precipitation₀` expects, so
# it is passed through unchanged. (The May-22 file used through tag v0.3.0, formerly v1.2.0, stored it linearly and was log-transformed here.)
#
# get_model draws from the global RNG (its land-water realisation), so every BRICK 2.0 model is built with
# `Random.seed!(BRICK20_SEED)` immediately before get_model. The hindcast drivers never read BRICK's own land water
# (it is zero before 2019 by calibration design); they add the observed series instead.

const BRICK20_FILE = ("brick20", "parameters_subsample_sneasybrick.csv")
const BRICK20_LABEL = "MimiBRICK v2.0.0, Zenodo 20592337 joint SNEASY-BRICK posterior"
const BRICK20_SEED = 2026

"""The BRICK 2.0 posterior subsample (10,000 rows). Errors with the fetch instruction if it is absent, and if
`antarctic_precip0` is positive in every row: that is a linear-units file (the May-22 one), not log space."""
function brick20_posterior()
    p = datapath(BRICK20_FILE...)
    isfile(p) || error("BRICK 2.0 posterior not found at $p; run tools/fetch_brick20_posterior.sh")
    df = CSV.read(p, DataFrame)
    all(>(0), df.antarctic_precip0) && error("$p: antarctic_precip0 is positive in every row, so it looks linear; " *
                                             "this arm expects log space (see the header of brick20.jl)")
    return df
end

"""
    brick20_model(; years=1850:2300, build_scenario=BUILD_SCENARIO) -> Mimi.Model

Stock MimiBRICK v2.0.0, built with the global RNG seeded to BRICK20_SEED immediately before get_model.
"""
function brick20_model(; years=1850:2300, build_scenario::AbstractString=BUILD_SCENARIO)
    Random.seed!(BRICK20_SEED)
    return MimiBRICK.get_model(ssprcp_scenario = build_scenario, start_year = first(years), end_year = last(years))
end

"""Set the prescribed GMST (K) and ocean heat content (1e22 J) of a BRICK 2.0 model."""
function set_brick20_forcing!(m, gmst, ohc)
    update_param!(m, :model_global_surface_temperature, gmst)
    update_param!(m, :thermal_expansion, :ocean_heat_interior, ohc)
    return m
end

"""Apply one BRICK 2.0 posterior row (MimiBRICK's own parameter names; precip already in log space)."""
function apply_brick20!(m, r)
    update_param!(m, :antarctic_ocean, :anto_α, r.anto_alpha)
    update_param!(m, :antarctic_ocean, :anto_β, r.anto_beta)
    for (col, sym) in ((:antarctic_s0, :ais_sea_level₀), (:antarctic_bed_height0, :ais_bedheight₀),
                       (:antarctic_slope, :ais_slope), (:antarctic_mu, :ais_μ),
                       (:antarctic_runoff_height0, :ais_runoffline_snowheight₀), (:antarctic_c, :ais_c),
                       (:antarctic_kappa, :ais_κ), (:antarctic_nu, :ais_ν), (:antarctic_flow0, :ais_iceflow₀),
                       (:antarctic_gamma, :ais_γ), (:antarctic_alpha, :ais_α),
                       (:antarctic_temp_threshold, :temperature_threshold), (:antarctic_lambda, :λ))
        update_param!(m, :antarctic_icesheet, sym, r[col])
    end
    update_param!(m, :antarctic_icesheet, :ais_precipitation₀, r.antarctic_precip0)
    update_param!(m, :glaciers_small_icecaps, :gsic_β₀, r.glaciers_beta0)
    update_param!(m, :glaciers_small_icecaps, :gsic_v₀, r.glaciers_v0)
    update_param!(m, :glaciers_small_icecaps, :gsic_s₀, r.glaciers_s0)
    update_param!(m, :glaciers_small_icecaps, :gsic_n,  r.glaciers_n)
    update_param!(m, :greenland_icesheet, :greenland_a,  r.greenland_a)
    update_param!(m, :greenland_icesheet, :greenland_b,  r.greenland_b)
    update_param!(m, :greenland_icesheet, :greenland_α,  r.greenland_alpha)
    update_param!(m, :greenland_icesheet, :greenland_β,  r.greenland_beta)
    update_param!(m, :greenland_icesheet, :greenland_v₀, r.greenland_v0)
    update_param!(m, :thermal_expansion, :te_α,  r.thermal_alpha)
    update_param!(m, :thermal_expansion, :te_s₀, r.thermal_s0)
    return m
end

"""BRICK 2.0 component variables: (component, variable)."""
const BRICK20_COMPONENTS = Dict(:ais => (:antarctic_icesheet, :ais_sea_level),
                                :glaciers => (:glaciers_small_icecaps, :gsic_sea_level),
                                :gis => (:greenland_icesheet, :greenland_sea_level),
                                :te => (:thermal_expansion, :te_sea_level),
                                :lws => (:landwater_storage, :lws_sea_level))

"""One BRICK 2.0 component of the last run, in cm relative to the mean over the `iref` indices."""
function brick20_sea_level(m, component::Symbol, iref)
    comp, var = BRICK20_COMPONENTS[component]
    v = m[comp, var]
    return 100 .* (v .- sum(v[iref]) / length(iref))
end
