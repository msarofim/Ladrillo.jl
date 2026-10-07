# The Ladrillo model: MimiBRICK v2.0.0 with the three-reservoir glacier component and the two-basin Greenland component
# in their slots, the Antarctic and thermal-expansion components stock, driven by prescribed GMST and ocean heat content.

const REF_PERIOD = (1995, 2014)          # reporting baseline for sea level (AR6-comparable)
const DRIVER_BASE = (1850, 1900)         # baseline of every temperature driver
const ANCHOR_YEARS = 11                  # observed years the anchor-preserving splice matches (the last 11)
const GIS_V0_M = 7.42                    # Greenland volume, m SLE (structural)
const GIS_G = 0.0                        # fraction of the 1850 commitment already realised (fixed)
"""T_ant at GMST = 0 in the Antarctic temperature map, preserved while the GMST->T_ant amplification is sampled."""
const AIS_TANT0 = -15.42 / 0.8365
"""Amplification-law averaging window (yr): the shape S is evaluated on the centred running mean of GMST."""
const GIS_SHAPE_WINDOW = 30
"""MimiBRICK build scenario. get_model's scenario only selects its bundled default forcing, which `LadrilloModel`
overwrites; it is fixed so any forcing can be run (MimiBRICK ships no file for most scenarios)."""
const BUILD_SCENARIO = "ssp245"
"""Mean GMST used to centre the sampled Antarctic precipitation coordinate (the calibrator's TBAR_ANT): AIS_TANT0 plus
the prior-mean amplification (1.09) times the 1900-2026 mean of the calibration forcing. A convention, not a fit."""
function tbar_ant()
    get!(_CACHE, :tbar_ant) do
        d = CSV.read(datapath("forcing", "fair_mean_gmst_ssp245harm.csv"), DataFrame)
        (-15.42 / 0.8365) + 1.09 * mean(d.gmst_C[(d.year .>= 1900) .& (d.year .<= 2026)])
    end
end

const _GLAC = :glaciers_small_icecaps
const _GIS  = :greenland_icesheet
const _AIS  = :antarctic_icesheet

"""Posterior columns assigned 1:1 to a model parameter: (column, component, parameter)."""
const DIRECT_PARAMS = [
    ("ais_ocean_temperature₀",   _AIS, :ais_ocean_temperature₀),
    ("antarctic_alpha",          _AIS, :ais_α),
    ("antarctic_nu",             _AIS, :ais_ν),
    ("antarctic_temp_threshold", _AIS, :temperature_threshold),
    ("antarctic_lambda",         _AIS, :λ),
    ("antarctic_gamma",          _AIS, :ais_γ),
    ("antarctic_kappa",          _AIS, :ais_κ),
    ("ais_mu",                   _AIS, :ais_μ),
    ("ais_bedheight0",           _AIS, :ais_bedheight₀),
    ("ais_slope",                _AIS, :ais_slope),
    ("ais_iceflow0",             _AIS, :ais_iceflow₀),
    ("ais_c",                    _AIS, :ais_c),
    ("anto_alpha",     :antarctic_ocean,   :anto_α),
    ("anto_beta",      :antarctic_ocean,   :anto_β),
    ("thermal_alpha",  :thermal_expansion, :te_α),
    ("gis_c1",      _GIS, :gis_c1),
    ("gis_c0",      _GIS, :gis_c0),
    ("gis_f",       _GIS, :gis_f),
    ("gis_alpha_f", _GIS, :gis_alpha_f),
    ("gis_beta_f",  _GIS, :gis_beta_f),
    ("gis_alpha_s", _GIS, :gis_alpha_s),
    ("gis_beta_s",  _GIS, :gis_beta_s),
]

"""
A built Ladrillo model on one forcing, plus the draw-independent pieces of the temperature drivers.
Construct with `LadrilloModel(gmst, ohc)` or `LadrilloModel(scenario)`.
"""
struct LadrilloModel
    m::Mimi.Model
    years::Vector{Int}
    iref::Vector{Int}                        # indices of REF_PERIOD
    gmst::Vector{Float64}                    # GMST as fed, K rel. 1850-1900
    gmst_rb::Vector{Float64}                 # GMST explicitly re-based to DRIVER_BASE
    glac_obs::Dict{String,Vector{Float64}}   # observed per-reservoir driver, padded to `years`
    glac_anchor::Dict{String,Float64}        # its mean over the anchor window
    glac_mask::BitVector                     # years with an observed glacier driver
    gmst_anchor::Float64                     # mean gmst_rb over the glacier anchor window
    gis_obs::Vector{Float64}                 # observed regional Greenland driver, padded
    gis_anchor::Float64
    gis_mask::BitVector
    gis_shape::Vector{Float64}               # amplification shape S(warming level) per year
    gis_shape_anchor::Float64                # mean(S * gmst_rb) over the Greenland anchor window
    tap::Bool
end

"""Greenland amplification shape S(dT): linear interpolation of the CMIP6-derived table, held flat outside it.
S(anchor dT) == 1 is asserted at first use: the law meets the calibration where the calibration was made."""
function gis_shape(dt::Real)
    tbl = get!(_CACHE, :gis_shape) do
        t = CSV.read(datapath("constants", "gis_amp_shape.csv"), DataFrame)
        meta = CSV.read(datapath("constants", "gis_amp_shape_meta.csv"), DataFrame)[1, :]
        (x = t.dt, y = t.S, anchor = Float64(meta.anchor_dt))
    end
    x, y = tbl.x, tbl.y
    d = clamp(Float64(dt), first(x), last(x))
    i = searchsortedlast(x, d)
    i >= length(x) && return Float64(y[end])
    w = (d - x[i]) / (x[i+1] - x[i])
    return Float64(y[i]) * (1 - w) + Float64(y[i+1]) * w
end

function _check_gis_shape_anchor()
    gis_shape(0.0)                                   # fills the cache
    a = _CACHE[:gis_shape].anchor
    abs(gis_shape(a) - 1.0) < 1e-9 || error("gis_amp_shape: S(anchor dT = $a) = $(gis_shape(a)) != 1")
end

"""Centred running mean over `w` years, shrinking at the ends."""
function running_mean(v::AbstractVector{<:Real}, w::Int)
    w <= 1 && return Float64.(v)
    n, lo, hi = length(v), (w - 1) ÷ 2, w ÷ 2
    return [mean(@view v[max(1, i - lo):min(n, i + hi)]) for i in 1:n]
end

_yearmap(path, col) = (d = CSV.read(path, DataFrame); Dict(Int(d[i, "year"]) => Float64(d[i, col]) for i in 1:nrow(d)))

"""Fixed per-reservoir nu (not sampled; the hindcast cannot identify it)."""
function glacier_nu()
    get!(_CACHE, :glacier_nu) do
        bc = CSV.read(datapath("constants", "glacier_block_constants.csv"), DataFrame)
        Dict(b => Float64(bc[findfirst(==(b), bc.block), "nu_anch_obsfit"]) for b in GLACIER_BLOCKS)
    end
end

"""The two parameters the posterior does not sample and no draw overwrites (Antarctic and thermal-expansion initial
sea level), from the calibration's central row."""
function fixed_params()
    get!(_CACHE, :fixed_params) do
        r = CSV.read(datapath("constants", "central_row.csv"), DataFrame)[1, :]
        (ais_sea_level0 = Float64(r.antarctic_s0), te_s0 = Float64(r.thermal_s0))
    end
end

"""
    LadrilloModel(gmst, ohc; years=1850:2300, lws=:observed, tap=true)

Build Ladrillo on prescribed GMST (K rel. 1850-1900) and ocean heat content (1e22 J), one value per year.
`lws` is :observed (projections) or :central (hindcast-side runs); `tap` switches the Greenland above-threshold
discharge channel on at the shipped TAP_CELL.
"""
function LadrilloModel(gmst::AbstractVector{<:Real}, ohc::AbstractVector{<:Real};
                       years=1850:2300, lws::Symbol=:observed, tap::Bool=true)
    years = collect(Int, years)
    y0, y1 = first(years), last(years)
    length(gmst) == length(years) || error("gmst has $(length(gmst)) values, need $(length(years))")
    length(ohc) == length(years) || error("ohc has $(length(ohc)) values, need $(length(years))")
    yi(y) = findfirst(==(y), years)
    _check_gis_shape_anchor()

    ibase = findall(y -> DRIVER_BASE[1] <= y <= DRIVER_BASE[2], years)
    gmst_rb = Float64.(gmst) .- mean(Float64.(gmst)[ibase])

    tgb = CSV.read(datapath("drivers", "glacier_block_temperature.csv"), DataFrame)
    last_obs = Int(maximum(tgb.year))
    last_obs <= y1 || error("observed glacier drivers end $last_obs, past the horizon $y1")
    anchor = (last_obs - ANCHOR_YEARS + 1):last_obs
    glac_obs, glac_anchor = Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    for b in GLACIER_BLOCKS
        d = Dict(Int(tgb[i, :year]) => Float64(tgb[i, b]) for i in 1:nrow(tgb))
        glac_obs[b] = [get(d, y, 0.0) for y in years]
        glac_anchor[b] = mean(d[y] for y in anchor)
    end

    tgz = CSV.read(datapath("drivers", "greenland_zone_temperature.csv"), DataFrame)
    gd = Dict(Int(tgz[i, :year]) => Float64(tgz[i, GIS_ZONE]) for i in 1:nrow(tgz))
    gis_last = Int(maximum(tgz.year))
    gis_last <= y1 || error("Greenland driver ends $gis_last, past $y1")
    ganch = (gis_last - ANCHOR_YEARS + 1):gis_last
    gis_obs = [get(gd, y, 0.0) for y in years]
    gis_anchor = mean(gd[y] for y in ganch)
    shape = gis_shape.(running_mean(gmst_rb, GIS_SHAPE_WINDOW))
    ianch = [yi(y) for y in ganch]
    shape_anchor = mean(shape[ianch] .* gmst_rb[ianch])

    m = MimiBRICK.get_model(ssprcp_scenario=BUILD_SCENARIO, start_year=y0, end_year=y1)
    replace!(m, _GLAC => ladrillo_glaciers)
    set_lws!(m, lws)
    replace!(m, _GIS => ladrillo_greenland)

    fp = fixed_params()
    update_param!(m, _AIS, :ais_sea_level₀, fp.ais_sea_level0)
    update_param!(m, :thermal_expansion, :te_s₀, fp.te_s0)
    nu = glacier_nu()
    for b in GLACIER_BLOCKS
        update_param!(m, _GLAC, Symbol("gic_nu_$b"), nu[b])
    end
    update_param!(m, _GLAC, :gic_sl0, 0.0)
    update_param!(m, :model_global_surface_temperature, Float64.(gmst))
    update_param!(m, :thermal_expansion, :ocean_heat_interior, Float64.(ohc))
    update_param!(m, _GIS, :gis_v0, GIS_V0_M)
    update_param!(m, _GIS, :gis_g, GIS_G)
    update_param!(m, _GIS, :gis_k_active, GIS_K_ACTIVE)
    update_param!(m, _GIS, :gis_k_high, GIS_K_HIGH)
    update_param!(m, _GIS, :gis_tap_onset, TAP_CELL.onset_K)
    update_param!(m, _GIS, :gis_tap_tau, TAP_CELL.tau_yr)
    update_param!(m, _GIS, :gis_tap_ramp_w, TAP_CELL.ramp_w_K)
    if tap
        maximum(gmst) > minimum(gmst) || error("tap on but the GMST series is constant")
        update_param!(m, _GIS, :gis_tap_v, TAP_CELL.V_m)
        update_param!(m, _GIS, :gis_tap_gmt, Float64.(gmst))
    else
        update_param!(m, _GIS, :gis_tap_v, 0.0)
        update_param!(m, _GIS, :gis_tap_gmt, zeros(length(years)))
    end

    return LadrilloModel(m, years, [yi(y) for y in REF_PERIOD[1]:REF_PERIOD[2]],
                         Float64.(gmst), gmst_rb, glac_obs, glac_anchor, years .<= last_obs,
                         mean(gmst_rb[[yi(y) for y in anchor]]),
                         gis_obs, gis_anchor, years .<= gis_last, shape, shape_anchor, tap)
end

"""Per-reservoir glacier driver: observations where they exist, then `amp * GMST` spliced on with the offset that
preserves the observed mean over the anchor window."""
glacier_driver(lm::LadrilloModel, block::AbstractString, amp::Real) =
    ifelse.(lm.glac_mask, lm.glac_obs[block], amp .* lm.gmst_rb .+ (lm.glac_anchor[block] - amp * lm.gmst_anchor))

"""Regional Greenland driver at amplification `amp`: observations, then `amp * S(warming level) * GMST`, anchor-
preserving. S = 1 at the calibration's own warming level."""
greenland_driver(lm::LadrilloModel, amp::Real) =
    ifelse.(lm.gis_mask, lm.gis_obs, amp .* lm.gis_shape .* lm.gmst_rb .+ (lm.gis_anchor - amp * lm.gis_shape_anchor))

"""
    apply_draw!(lm, row)

Apply one posterior draw (a row of `posterior()`) to the model without running it: the direct parameters, the derived
Antarctic quantities (log P0 from the sampled coordinate; the runoff line from its onset temperature; the temperature
map with T_ant(GMST = 0) preserved), the per-reservoir glacier parameters and drivers, the Greenland driver and the
high-basin rate scale.
"""
function apply_draw!(lm::LadrilloModel, row)
    m = lm.m
    for (col, comp, sym) in DIRECT_PARAMS
        update_param!(m, comp, sym, Float64(row[col]))
    end
    update_param!(m, _AIS, :ais_precipitation₀,
                  Float64(row["ais_precip_u"]) - Float64(row["antarctic_kappa"]) * tbar_ant())
    update_param!(m, _GIS, :greenland_surface_temperature, Float64.(greenland_driver(lm, Float64(row["gis_amp"]))))
    update_param!(m, _GIS, :gis_s_high, 10.0^Float64(row["gis_s_high"]))
    for b in GLACIER_BLOCKS
        update_param!(m, _GLAC, Symbol("gic_a_$b"), Float64(row["gic_a_$b"]))
        update_param!(m, _GLAC, Symbol("gic_b_$b"), Float64(row["gic_b_$b"]))
        update_param!(m, _GLAC, Symbol("gic_T_off_$b"), Float64(row["gic_T_off_$b"]))
        update_param!(m, _GLAC, Symbol("gic_kappa_$b"), 10.0^Float64(row["gic_log10_kappa_$b"]))
        update_param!(m, _GLAC, Symbol("glacier_surface_temperature_$b"),
                      glacier_driver(lm, b, Float64(row["gic_amp_$b"])))
    end
    update_param!(m, _AIS, :ais_runoffline_snowheight₀, -Float64(row["ais_runoff_Ton"]) * Float64(row["ais_c"]))
    amp = Float64(row["ais_gmst_amp"])
    update_param!(m, _AIS, :ais_temperature_coefficient, 1.0 / amp)
    update_param!(m, _AIS, :ais_temperature_intercept, -AIS_TANT0 / amp)
    return lm
end

"""`apply_draw!` then run."""
run_draw!(lm::LadrilloModel, row) = (apply_draw!(lm, row); run(lm.m); lm)

"""Sea-level components readable with `sea_level`. `total` is MimiBRICK's summed global sea level
(glaciers + Greenland + Antarctica + thermal expansion + land water)."""
const COMPONENTS = Dict{Symbol,Tuple{Symbol,Symbol}}(
    :glaciers   => (_GLAC, :gsic_sea_level),
    :gsic_r19   => (_GLAC, :gsic_r19),
    :gsic_slowg => (_GLAC, :gsic_slowg),
    :gsic_fastg => (_GLAC, :gsic_fastg),
    :gsic_hind  => (_GLAC, :gsic_hind),
    :gis        => (_GIS,  :greenland_sea_level),
    :ais        => (_AIS,  :ais_sea_level),
    :te         => (:thermal_expansion, :te_sea_level),
    :lws        => (:landwater_storage, :lws_sea_level),
    :total      => (:global_sea_level,  :sea_level_rise),
)

"""
    sea_level(lm, component) -> Vector{Float64}

One component of the last run, in cm relative to REF_PERIOD.
"""
function sea_level(lm::LadrilloModel, component::Symbol)
    haskey(COMPONENTS, component) || error("unknown component :$component; have $(sort(collect(keys(COMPONENTS))))")
    comp, var = COMPONENTS[component]
    x = [ismissing(e) ? NaN : 100.0 * Float64(e) for e in lm.m[comp, var]]
    return x .- mean(x[lm.iref])
end
