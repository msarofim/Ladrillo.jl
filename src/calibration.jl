# The calibration (L27): the log-posterior Ladrillo's posterior was sampled from, its parameter table, the start points
# and the proposal-covariance seed. scripts/calibrate.jl samples it.
#
# A port of the research calibrator (SLR-RFF-BRICK julia/calibrate_mcmc_ext.jl, run with the flags of run_mcmc_L27.sh,
# commit 43d84aa) reduced to that one configuration:
#   --gis-ordered --gis-basins2 --amp-mu=1.09 --amp-sigma=0.180 --paleo-priors --no-delta --no-d2-gsic
#   --obs-corr-len=100 --toff-lo=-4 --precip-reparam --cut-fastdyn --fix-gamma --no-ledger
# The model is the package's own LadrilloModel on 1850-2026 (calibration forcing ssp245harm, land water :central,
# Greenland tap off). Every likelihood and prior term keeps the research code's arithmetic and order of accumulation, so
# on the reference runtime the port reproduces the production chain bit for bit (test/calibration/).
#
# PARAMETERS (42 physical + 8 noise, in sampling order; `parameter_table()`):
#   Antarctica   ais_ocean_temperature0, antarctic_{alpha,nu,kappa} and anto_{alpha,beta} (paleo-ensemble marginals),
#                the GMST->Antarctic amplification ais_gmst_amp ~ N(1.09, 0.18), and the seven DAIS geometry parameters
#                under their joint paleo prior. Two of those are sampled in identified coordinates: the runoff onset
#                ais_runoff_Ton = -h0/c, and ais_precip_u = log P0 + kappa * TBAR_ANT (log P0 is derived per draw and
#                the paleo prior and bounds apply to it). lambda, T_crit and gamma are not sampled (HELD, below).
#   Greenland    gis_c1, gis_c0, gis_f, gis_alpha_f, gis_beta_f; the slow channel as (gis_slow_ell, gis_slow_w) at the
#                driver's 2015-2024 mean; gis_s_high (log10 rate scale of the high basin); gis_amp (prior only).
#   Glaciers     per reservoir (R19, SLOWG, FASTG): gic_a, gic_b, gic_T_off, gic_log10_kappa (prior centred on a
#                log-linear function of the reservoir's amplification), gic_amp (prior only); gic_u_unch (uncharted ice).
#   Thermal exp. thermal_alpha (flat on [0.05, 0.30]); d2_steric_{1,2}, a two-term discrepancy orthogonal to the
#                steric driver shape.
#   Noise        sd and rho of an AR(1) residual for each scored series (ais, gsic, gis, steric).
#
# LIKELIHOOD. The four component series 1900-2025/26 (each an AR(1) residual plus the observational band as an error
# correlated over OBS_CORR_LEN years); Antarctic surface mass balance 1979-2008; the glacier inventory in 2000; the
# 1850-1900 glacier loss (with the two pre-1901 set-asides marginalised over their priors); the per-reservoir
# GlacierMIP3 committed-loss rungs; the GlaMBIE 2000-2024 SLOWG/FASTG partition and R19 rate; the Mouginot 2019
# surface share of Greenland's extra loss; and the Mouginot active-basin share in 2002-2011 and 2012-2018. The total
# (GMSL) is not scored.
#
# THE DRIVERS ARE BUILT ONCE, at the prior centres of the amplifications. The amplification-dependent part of every
# driver is the splice after the observed record (2025-2026 here), which no glacier term reads; the Greenland target's
# 2025 value is the one year that sees it. So gic_amp and gis_amp are carried by their priors into the projections,
# which rebuild the drivers per draw (model.jl), and are not estimated.

module Calibration

using CSV, DataFrames, Mimi, Statistics
import Distributions: Normal, MvNormal, logpdf, truncated
import LinearAlgebra: Symmetric, Diagonal, I, dot, norm, isposdef, diag
import Random
import RobustAdaptiveMetropolisSampler: RAM_sample
using ..Ladrillo: LadrilloModel, datapath, mean_forcing, funch_unit, tbar_ant, GLACIER_BLOCKS, FIT_REF, AIS_TANT0

const TAG = "L27"
const Y0, Y1 = 1850, 2026
const FORCING = "ssp245harm"
const BLOCKS = collect(GLACIER_BLOCKS)
const SERIES = (:ais, :gsic, :gis, :steric)              # scored component series; the total is not scored
const FIT_START = 1900

## AR(1) noise and the observational band
const OBS_CORR_LEN = 100.0      # yr, e-folding length of the band's error correlation (--obs-corr-len=100)
const RHO_MAX = 0.99            # hard upper bound on every AR(1) rho
const NOISE_SD_PRIOR = 5.0      # cm, half-normal prior on every AR(1) sd
const BAND_Z, BAND_FLOOR = 1.645, 0.05   # per-year sigma = half the 90% band / 1.645, floored at 0.05 cm

## Glaciers
"""Per-reservoir amplification prior (mean, sd, lo, hi): glacier-area warming per K of GMST."""
const AMP_PRIOR = Dict("R19" => (0.72, 0.15, 0.58, 0.88),
                       "SLOWG" => (2.50, 0.45, 1.80, 3.50),
                       "FASTG" => (1.45, 0.15, 1.33, 1.82))
const K10_SIG = 0.114           # sd of the log10 kappa prior (+-30 % at 1 sigma)
const FIT_BASIS = "obsfit"      # glacier_block_constants.csv column suffix for the b, T_off centres and nu
const TOFF_LO, TOFF_HI = -4.0, 1.0       # --toff-lo=-4
const GMIP_LEVELS = [1.2, 1.5, 2.0, 3.0] # GlacierMIP3 warming levels, K
const RUNG_CORR = 0.6
const RUNG_D2_N8 = 2.847        # expected range of 8 normal draws in sd (d2(8)): the rung sd is the 8-model range / d2
const RUNG_SIG_SCALE = 2.0 / RUNG_D2_N8
const GLAMBIE_WIN = (2000, 2024)
const R19_RATE_SD = 0.11615     # mm/yr, GlaMBIE R19 rate sd with the per-year errors fully serially correlated
const GLAMBIE_SHARE_SD = 0.05
const GLAMBIE_TOT_FLOOR = 1e-9
const INV_V_M, INV_SIGMA_M, INV_YEAR = 0.290, 0.060, 2000   # glacier inventory (RGI 1-18 minus 5, plus 19), m SLE
const M19_MU_M, M19_SIGMA_M = 0.020, 0.009                  # 1850-1900 glacier loss (Leclercq family), m SLE
const LEDGER_SHIFT_M = 0.015    # mean of the two pre-1901 set-asides (uniform[0, 25] mm and N(2.5, 2.0) mm), m
const LEDGER_SD_M = (0.0072, 0.0020)   # their sds, m, added in quadrature

## Greenland
const GIS_ZONE, GIS_AMP_WINDOW = "south", "full"
const GIS_TBAR_WIN = (2015, 2024)
const GIS_TBAR_REPARAM_ANCHOR = 1.963    # K, the driver mean the (ell, w) coordinates were chosen at
const GIS_NATIVE_MU = (alpha_s = 0.0070727, beta_s = 0.0010)   # the slow channel's prior centres in native coordinates
const GIS_ELL_SD = 1.0
const MOUG_SHARE, MOUG_SHARE_SD = 0.735, 0.05
const MOUG_REF_WIN, MOUG_LATE_WIN = (1972, 1990), (2000, 2018)
const GISB_WINS = ((2002, 2011), (2012, 2018))
"""Mouginot 2019 sector shares of the mean loss rate in GISB_WINS. The active basin is south + mid; high follows by
sum-to-one and is not scored."""
const GISB_SHARE3 = ((south = 0.592, mid = 0.207, high = 0.201),
                     (south = 0.554, mid = 0.262, high = 0.183))
const GISB_SHARE_ACTIVE = Tuple(sh.south + sh.mid for sh in GISB_SHARE3)
const GISB_SHARE_SD = 0.05
const GISB_TOT_FLOOR = 1e-12

## Antarctica
const AMP_MU, AMP_SIGMA = 1.09, 0.180
"""DAIS parameters with paleo-ensemble marginal priors (paleo_dais_marginals.csv: mean, sd, lo, hi)."""
const PALEO_PRIOR = Set(["anto_alpha", "anto_beta", "antarctic_alpha", "antarctic_nu", "antarctic_kappa"])
"""
Not sampled; held during the calibration at the CENTRAL ROW's values (data/constants/central_row.csv):
lambda 0.01365, T_crit -15.14, gamma 2.038.

THESE ARE NOT THE PALEO MEDIANS, although the research calibrator's --cut-fastdyn and --fix-gamma flags say so. It set
the paleo medians (lambda 0.00991, T_crit -15.67, gamma 2.834) and then loaded the central row over them
(calibrate_mcmc_ext.jl:1605-1607, then update_brick_nu3!(m, medoid, ...) at :1627), so the production L27 chain was
sampled with the central-row values. Found 2026-11-11 by this port: with the paleo medians the log-posterior at the four
start rows misses production by -0.015 to +0.034; with the central row it matches all four exactly.
- lambda and T_crit are inert in the hindcast (the fast-dynamics threshold is never crossed): either pair gives a
  bit-identical log-posterior.
- gamma is not inert. posterior() attaches the paleo median (2.834) for projections and the posterior predictive, so
  those products use a gamma the chain was not conditioned on. Measured on 500 draws: hindcast AIS moves <= 0.015 cm,
  projections <= 0.13 cm through 2100 and up to 3.2 cm in the 2300 AIS median (SSP5-8.5).
"""
const HELD = [("antarctic_lambda", :λ), ("antarctic_temp_threshold", :temperature_threshold), ("antarctic_gamma", :ais_γ)]
const GEO_NAMES = ["ais_mu", "ais_bedheight0", "ais_slope", "ais_iceflow0", "ais_precip0_LOG", "ais_runoff_Ton", "ais_c"]
const GEO_SYMS = [:ais_μ, :ais_bedheight₀, :ais_slope, :ais_iceflow₀, :ais_precipitation₀, :ais_runoffline_snowheight₀,
                  :ais_c]
const SMB_WIN = (1979, 2008)
const SMB_TARGET_GT = 2098.0 * (10.92 / 12.295)   # Rignot 2019 grounded SMB, area-scaled to the DAIS disc
const SMB_SIGMA_GT = 133.0 * (10.92 / 12.295)
const M3ICE_TO_GT = 917.0 / 1e12

## Sampler set-up
const D2_BASIS_N, D2_BASIS_SD = 2, 0.5
const SEEDS_CANON = [2026, 2027, 2028, 2029]       # chain i starts from row i of the starts file
const SEEDS_REPLICATE = [3026, 3027, 3028, 3029]   # the same start rows on another RNG stream
const STARTS_FILE = "overdispersed_starts_L27.csv"
const ADCOV_FILE = "adapted_cov_L26_named.csv"
const PROP_SCALE, GEO_PROP_SCALE = 0.1, 0.02
const NOISE_PROP = (0.3, 0.1)
const OPT_ALPHA = 0.234
"""Floors on the proposal-seed sd of the AIS geometry block and the basin scale: a seed that hands a parameter its
neighbour's (far smaller) variance freezes it for the whole run while global acceptance looks healthy."""
const GEO_SEED_FLOOR = Dict("ais_c" => 0.05, "ais_mu" => 3e-3, "ais_bedheight0" => 0.03, "ais_slope" => 1e-8,
                            "ais_iceflow0" => 3e-4, "ais_precip0_LOG" => 1e-3, "ais_runoff_Ton" => 5e-4,
                            "gis_s_high" => 1e-3)

const _AIS, _GIS, _GLAC = :antarctic_icesheet, :greenland_icesheet, :glaciers_small_icecaps
const LIK = :likelihood_only

"""AR(1) residual with stationary variance sd^2/(1-rho^2) plus the observational band eps as an error correlated with
e-folding length L yr (Ruckert et al. 2017, the MimiBRICK form)."""
function hetero_logl_ar1(res::Vector{Float64}, σ::Float64, ρ::Float64, ϵ::Vector{Float64}, L::Float64)
    n = length(res)
    σp = σ^2/(1-ρ^2)
    H  = abs.(collect(1:n)' .- collect(1:n))
    Σ  = L > 0 ? σp .* ρ.^H .+ (ϵ * ϵ') .* exp.(-H ./ L) :
                 σp .* ρ.^H .+ Diagonal(ϵ.^2)
    return logpdf(MvNormal(Symmetric(Σ)), res)
end

"""Orthonormal (unit-RMS) discrepancy basis: powers of scaled time, Gram-Schmidt'd against `protect` (plain inner
product) and each other. (nyear x D2_BASIS_N)"""
function d2_basis(years, protect::Vector{Vector{Float64}})
    n = length(years)
    ip = (u, v) -> dot(u, v)
    x = 2 .* (Float64.(years) .- minimum(years)) ./ (maximum(years) - minimum(years)) .- 1
    cols = Vector{Vector{Float64}}()
    base = Vector{Vector{Float64}}()
    for u0 in protect
        u = copy(u0)
        for w in base
            d = ip(w, w); d > 1e-12 && (u = u .- (ip(u, w) / d) .* w)
        end
        norm(u) > 1e-10 && push!(base, u)
    end
    for k in 1:(D2_BASIS_N + length(protect) + 2)
        length(cols) == D2_BASIS_N && break
        v = x .^ k
        for u in base
            d = ip(u, u)
            d > 1e-12 && (v = v .- (ip(v, u) / d) .* u)
        end
        rms = sqrt(sum(v .^ 2) / n)
        rms < 1e-8 && continue
        v = v ./ rms
        push!(cols, v); push!(base, v)
    end
    length(cols) == D2_BASIS_N || error("d2_basis: only $(length(cols)) of $D2_BASIS_N independent columns")
    return hcat(cols...)
end

"""One scored series: fit years, model indices, observations, per-year band sigma."""
const Series = NamedTuple{(:years, :myi, :obs, :ϵ), Tuple{Vector{Int}, Vector{Int}, Vector{Float64}, Vector{Float64}}}

"""
    Problem

The calibration's fixed pieces: the built model, the parameter table (sampling order), the scored series, the
discrepancy basis and every constant the log-posterior reads. Build with `problem()`.
"""
struct Problem
    lm::LadrilloModel
    years::Vector{Int}
    ib::Vector{Int}                         # FIT_REF indices (the re-reference window)
    names::Vector{String}                   # the 42 physical parameters, then sd_/rho_ per series
    comp::Vector{Symbol}
    sym::Vector{Symbol}
    mu::Vector{Float64}
    sigma::Vector{Float64}
    lo::Vector{Float64}
    hi::Vector{Float64}
    np::Int
    nk::Int
    set_direct::Vector{Int}                 # set 1:1 on a model parameter
    prior_gauss::Vector{Int}                # independent N(mu, sigma) prior
    ix::Dict{String,Int}
    geo_idx::UnitRange{Int}
    geo_mu::Vector{Float64}
    geo_sd::Vector{Float64}
    geo_prior::MvNormal
    precip_lo::Float64                      # paleo bounds on the derived log P0
    precip_hi::Float64
    tbar_ant::Float64
    gis_tbar::Float64
    series::NamedTuple{SERIES, NTuple{4, Series}}
    d2_steric::Matrix{Float64}
    funch::Vector{Float64}
    k10_pts::Dict{String,NTuple{2,NTuple{2,Float64}}}
    s2020::Dict{String,Float64}
    rung_y::Dict{String,Vector{Float64}}
    rung_ci::Dict{String,Matrix{Float64}}
    r19_rate_mu::Float64
    glambie_fast_share::Float64
    held::Dict{String,Float64}              # lambda, T_crit, gamma during the calibration (HELD)
    gis_amp::Float64                        # the Greenland amplification the driver is built at (prior mean)
    i::NamedTuple                           # year indices used by point terms
end

idx(p::Problem, name::AbstractString) = p.ix[name]
"""Parameter names in sampling order (the chain's columns before log_post)."""
param_names(p::Problem) = p.names

"""κ-prior centre (log10) as a function of the reservoir amplification: log-linear between the two anchored solves."""
function k10c(p::Problem, b, amp)
    (a1, k1), (a2, k2) = p.k10_pts[b]
    return k1 + (k2 - k1) * (amp - a1) / (a2 - a1)
end

"""log P0 from the sampled u = log P0 + kappa * TBAR_ANT."""
precip_log(p::Problem, θ) = θ[p.ix["ais_precip_u"]] - θ[p.ix["antarctic_kappa"]] * p.tbar_ant

"""
    problem() -> Problem

Build the calibration: forcing, targets, model, parameter table, discrepancy basis and point-term constants. Asserts
the identities the reparameterisations rest on.
"""
function problem()
    years = collect(Y0:Y1)
    yidx(y) = findfirst(==(y), years)
    ib = [yidx(y) for y in FIT_REF[1]:FIT_REF[2]]
    gmst, ohc = mean_forcing(FORCING; years=years)
    tg = CSV.read(datapath("observations", "calibration_targets.csv"), DataFrame)
    bcdf = CSV.read(datapath("constants", "glacier_block_constants.csv"), DataFrame)
    bcrow(b) = bcdf[findfirst(==(b), bcdf.block), :]

    ## ---- glacier drivers: observed reservoir temperature, then amp x GMST (anchor-preserving), amp at the prior centre
    amp_b = Dict(b => AMP_PRIOR[b][1] for b in BLOCKS)
    k10_pts = Dict(b => ((Float64(bcrow(b).amp_regchar), log10(Float64(bcrow(b).kappa_anch_regchar))),
                         (Float64(bcrow(b).amp_obsfit), log10(Float64(bcrow(b).kappa_anch_obsfit)))) for b in BLOCKS)
    k10c0(b, amp) = (((a1, k1), (a2, k2)) = k10_pts[b]; k1 + (k2 - k1) * (amp - a1) / (a2 - a1))
    kap_anch = Dict(b => 10.0^k10c0(b, amp_b[b]) for b in BLOCKS)
    tgb = CSV.read(datapath("drivers", "glacier_block_temperature.csv"), DataFrame)
    tgb_last = Int(maximum(tgb.year))
    gmst_rb = gmst .- mean(gmst[findall(y -> 1850 <= y <= 1900, years)])
    tg3 = let gmd = Dict(zip(years, gmst_rb)), anchor = (tgb_last-10):tgb_last
        d = Dict{String,Vector{Float64}}()
        for b in BLOCKS
            obsd = Dict(Int(tgb[i, :year]) => Float64(tgb[i, b]) for i in 1:nrow(tgb))
            off = mean(obsd[y] for y in anchor) - amp_b[b] * mean(gmd[y] for y in anchor)
            d[b] = [y <= tgb_last ? obsd[y] : amp_b[b]*gmd[y] + off for y in years]
        end
        d
    end

    ## ---- Greenland: the amp prior, the (ell, w) anchor and the regional driver at the prior mean
    gap = let df = CSV.read(datapath("calibration", "gis_amp_prior.csv"), DataFrame)
        i = findfirst(r -> r.zone == GIS_ZONE && r.window == GIS_AMP_WINDOW, eachrow(df))
        isnothing(i) && error("gis_amp_prior.csv has no row zone=$GIS_ZONE window=$GIS_AMP_WINDOW")
        r = df[i, :]
        (mean=Float64(r.mean), sd=Float64(r.sd), lo=Float64(r.lo), hi=Float64(r.hi))
    end
    tgz = CSV.read(datapath("drivers", "greenland_zone_temperature.csv"), DataFrame)
    tgz_last = Int(maximum(tgz.year))
    gis_tbar = mean(Float64(tgz[i, GIS_ZONE]) for i in 1:nrow(tgz)
                    if GIS_TBAR_WIN[1] <= Int(tgz[i, :year]) <= GIS_TBAR_WIN[2])
    abs(gis_tbar - GIS_TBAR_REPARAM_ANCHOR) < 5e-3 ||
        error("Greenland driver mean $gis_tbar over $GIS_TBAR_WIN is not the $GIS_TBAR_REPARAM_ANCHOR K the (ell, w) " *
              "coordinates were chosen at; greenland_zone_temperature.csv has changed")
    ell_mu = log(GIS_NATIVE_MU.alpha_s * gis_tbar + GIS_NATIVE_MU.beta_s)
    w_mu = GIS_NATIVE_MU.alpha_s * gis_tbar / (GIS_NATIVE_MU.alpha_s * gis_tbar + GIS_NATIVE_MU.beta_s)
    let aM = w_mu * exp(ell_mu) / gis_tbar, bM = (1 - w_mu) * exp(ell_mu)    # the coordinates are a pure change
        (abs(aM - GIS_NATIVE_MU.alpha_s) < 1e-12 && abs(bM - GIS_NATIVE_MU.beta_s) < 1e-12) ||
            error("(ell, w) prior centres map to (alpha_s $aM, beta_s $bM), not $(GIS_NATIVE_MU)")
    end
    tgis = let obsd = Dict(Int(tgz[i, :year]) => Float64(tgz[i, GIS_ZONE]) for i in 1:nrow(tgz)),
               gmd = Dict(zip(years, gmst_rb)), anchor = (tgz_last-10):tgz_last
        off = mean(obsd[y] for y in anchor) - gap.mean * mean(gmd[y] for y in anchor)
        [y <= tgz_last ? obsd[y] : gap.mean * gmd[y] + off for y in years]
    end

    ## ---- scored series
    rowof(y) = findfirst(==(y), tg.year)
    function make_series(col, lo, hi)
        fy = Int[]
        for i in 1:nrow(tg)
            v = tg[i, col]
            (tg.year[i] >= FIT_START && !ismissing(v) && !isnan(Float64(v))) && push!(fy, Int(tg.year[i]))
        end
        sort!(fy)
        fy == collect(fy[1]:fy[end]) || error("series $col has a year gap (the AR(1) assumes unit spacing)")
        ri = [rowof(y) for y in fy]
        ϵ = max.((Float64.(tg[ri, hi]) .- Float64.(tg[ri, lo])) ./ (2*BAND_Z), BAND_FLOOR)
        return (years=fy, myi=[yidx(y) for y in fy], obs=Float64.(tg[ri, col]), ϵ=ϵ)
    end
    S = (ais = make_series(:ais, :ais_lo, :ais_hi), gsic = make_series(:gsic, :gsic_lo, :gsic_hi),
         gis = make_series(:gis, :gis_lo, :gis_hi), steric = make_series(:steric, :steric_lo, :steric_hi))
    ## the glacier target with the RGI-19 seam removed: equal to the raw target before 2019 (the model scope is SLOWG+FASTG)
    gadj = CSV.read(datapath("observations", "calibration_targets_gsic_adj.csv"), DataFrame)
    Int.(gadj.year) == S.gsic.years || error("calibration_targets_gsic_adj.csv year grid differs from the glacier target")
    let pre = S.gsic.years .< 2019
        maximum(abs.(Float64.(gadj.gsic_adj)[pre] .- S.gsic.obs[pre])) < 1e-6 ||
            error("the seam-adjusted glacier target differs from the raw target before 2019")
    end
    S.gsic.obs .= Float64.(gadj.gsic_adj)

    ## ---- the steric discrepancy basis, orthogonal to the constant and to the steric driver shape S(t) = OHC(t) + c
    te_shape = Float64.(ohc)[S.steric.myi]
    B = d2_basis(S.steric.years, [ones(length(S.steric.years)), copy(te_shape)])
    for k in 1:D2_BASIS_N
        v = B[:, k]
        abs(sum(v) / length(v)) < 1e-9 || error("steric discrepancy column $k is not mean-zero")
        abs(sqrt(sum(v .^ 2) / length(v)) - 1) < 1e-8 || error("steric discrepancy column $k is not unit-RMS")
        abs(dot(v, te_shape)) / (norm(v) * norm(te_shape)) < 1e-8 ||
            error("steric discrepancy column $k is not orthogonal to the steric shape (degenerate with thermal_alpha)")
    end
    abs(dot(B[:, 1], B[:, 2])) / (norm(B[:, 1]) * norm(B[:, 2])) < 1e-8 || error("steric discrepancy columns not orthogonal")

    ## ---- the parameter table, in sampling order
    pal = CSV.read(datapath("constants", "paleo_dais_marginals.csv"), DataFrame)
    palrow(n) = pal[findfirst(==(n), pal.param), :]
    held = let cr = CSV.read(datapath("constants", "central_row.csv"), DataFrame)[1, :]
        Dict(n => Float64(cr[n]) for (n, _) in HELD)
    end
    T = NamedTuple{(:name, :comp, :sym, :μ, :σ, :lo, :hi), Tuple{String, Symbol, Symbol, Float64, Float64, Float64, Float64}}
    F = T[]
    add!(n, c, s, μ, σ, lo, hi) = push!(F, (name=n, comp=c, sym=s, μ=Float64(μ), σ=Float64(σ), lo=Float64(lo), hi=Float64(hi)))
    paleo!(n, c, s) = (n in PALEO_PRIOR || error("$n has no paleo prior");
                       r = palrow(n); add!(n, c, s, r.mean, r.sd, r.lo, r.hi))
    add!("ais_ocean_temperature₀", _AIS, :ais_ocean_temperature₀, 0.72, 0.50, 0.50, 2.00)
    paleo!("antarctic_alpha", _AIS, :ais_α)
    paleo!("antarctic_nu", _AIS, :ais_ν)
    paleo!("anto_alpha", :antarctic_ocean, :anto_α)
    paleo!("anto_beta", :antarctic_ocean, :anto_β)
    ## Greenland: weak priors centred on the offline A+B fit at g = 0 (to the same target), level parameters in m SLE
    add!("gis_c1", _GIS, :gis_c1, 0.032766, 0.050, 0.0, 4.0)
    add!("gis_c0", _GIS, :gis_c0, 0.040429, 0.100, 0.0, 4.0)
    add!("gis_f", _GIS, :gis_f, 0.782569, 0.30, 0.02, 0.98)
    add!("gis_alpha_f", _GIS, :gis_alpha_f, 0.0028487, 0.020, 0.0, 0.5)
    add!("gis_beta_f", _GIS, :gis_beta_f, 0.0073684, 0.050, 1e-6, 0.5)
    add!("gis_slow_ell", LIK, :none, ell_mu, GIS_ELL_SD, ell_mu - 4*GIS_ELL_SD, ell_mu + 4*GIS_ELL_SD)
    add!("gis_slow_w", LIK, :none, w_mu, 1e3, 0.0, 1.0)
    add!("gis_amp", LIK, :none, gap.mean, gap.sd, gap.lo, gap.hi)
    add!("gis_s_high", _GIS, :gis_s_high, 0.0, 0.5, -2.0, 2.0)
    add!("thermal_alpha", :thermal_expansion, :te_α, 0.16, 1.0e3, 0.05, 0.30)   # MimiBRICK's Uniform(0.05, 0.3)
    for b in BLOCKS
        r = bcrow(b)
        a_lo = max(1.5*Float64(r.S2020_data), Float64(r.a0) - 3.5*Float64(r.a0_sig), 0.01)
        add!("gic_a_$b", _GLAC, Symbol("gic_a_$b"), r.a0, r.a0_sig, a_lo, Float64(r.a0) + 4.0*Float64(r.a0_sig))
        add!("gic_b_$b", _GLAC, Symbol("gic_b_$b"), r["b_fit_$(FIT_BASIS)"], 10.0, 0.05, 3.0)
        add!("gic_T_off_$b", _GLAC, Symbol("gic_T_off_$b"), r["T_off_fit_$(FIT_BASIS)"], 10.0, TOFF_LO, TOFF_HI)
        c1, c2 = k10c0(b, AMP_PRIOR[b][3]), k10c0(b, AMP_PRIOR[b][4])
        add!("gic_log10_kappa_$b", _GLAC, Symbol("gic_kappa_$b"), log10(kap_anch[b]), K10_SIG,
             min(c1, c2) - 1.0, max(c1, c2) + 1.0)
    end
    for b in BLOCKS
        μa, σa, loa, hia = AMP_PRIOR[b]
        add!("gic_amp_$b", LIK, :none, μa, σa, loa, hia)
    end
    add!("gic_u_unch", LIK, :none, 28.15, 1.0e3, 14.5, 41.8)
    for k in 1:D2_BASIS_N
        add!("d2_steric_$k", LIK, :none, 0.0, D2_BASIS_SD, -5*D2_BASIS_SD, 5*D2_BASIS_SD)
    end
    paleo!("antarctic_kappa", _AIS, :ais_κ)
    add!("ais_gmst_amp", _AIS, :ais_temperature_coefficient, AMP_MU, AMP_SIGMA, AMP_MU - 3*AMP_SIGMA, AMP_MU + 3*AMP_SIGMA)
    ## the DAIS geometry block, LAST: joint paleo prior on z = (theta - mean)/sd (paleo_geo_prior_ton.csv)
    gl = [split(strip(l), ',') for l in readlines(datapath("constants", "paleo_geo_prior_ton.csv"))
          if !startswith(l, "#") && !isempty(strip(l))]
    grow(tag) = [parse(Float64, x) for x in first(l for l in gl if l[1] == tag)[2:end]]
    String.(first(l for l in gl if l[1] == "names")[2:end]) == GEO_NAMES ||
        error("paleo_geo_prior_ton.csv columns are not $(GEO_NAMES)")
    geo_mu, geo_sd, glo, ghi = grow("mean"), grow("sd"), grow("lo"), grow("hi")
    geo_c = Symmetric(permutedims(reduce(hcat, [[parse(Float64, x) for x in l[2:end]] for l in gl if l[1] == "corr"])))
    for i in eachindex(GEO_SYMS)
        add!(GEO_NAMES[i], _AIS, GEO_SYMS[i], geo_mu[i], geo_sd[i], glo[i], ghi[i])
    end
    geo_idx = (length(F)-length(GEO_SYMS)+1):length(F)
    ## precipitation in its identified coordinate u = log P0 + kappa * TBAR_ANT; the bounds apply to the derived log P0
    tbar = tbar_ant()
    tbar == AIS_TANT0 + AMP_MU * mean(gmst[[yidx(y) for y in FIT_START:Y1]]) ||
        error("tbar_ant() is not AIS_TANT0 + AMP_MU x the $FIT_START-$Y1 mean of the calibration forcing")
    ip = geo_idx[findfirst(==("ais_precip0_LOG"), GEO_NAMES)]
    precip_lo, precip_hi = F[ip].lo, F[ip].hi
    let f = F[ip], κμ = F[findfirst(k -> k.name == "antarctic_kappa", F)].μ
        F[ip] = (name="ais_precip_u", comp=f.comp, sym=f.sym, μ=f.μ + κμ * tbar, σ=f.σ, lo=-1e9, hi=1e9)
    end

    np = length(F)
    names = vcat([f.name for f in F], vcat([["sd_$s", "rho_$s"] for s in SERIES]...))
    ix = Dict(n => i for (i, n) in enumerate(names))
    length(ix) == length(names) || error("duplicate parameter names")
    ## likelihood-only and derived parameters are never set 1:1 on the model
    derived = Set(vcat([ix["gic_log10_kappa_$b"] for b in BLOCKS], [ix["gic_amp_$b"] for b in BLOCKS],
                       [ix["d2_steric_$k"] for k in 1:D2_BASIS_N],
                       [ix[n] for n in ("gic_u_unch", "gis_slow_ell", "gis_slow_w", "gis_amp", "gis_s_high",
                                        "ais_gmst_amp", "ais_runoff_Ton", "ais_precip_u")]))
    set_direct = [k for k in 1:np if !(k in derived)]
    all(F[k].comp !== LIK && F[k].sym !== :none for k in set_direct) || error("a likelihood-only parameter is set on the model")
    ## the kappa priors are amplification-dependent (k10c) and the geometry block has the joint prior
    prior_gauss = [k for k in 1:np if !(k in geo_idx) && !(k in [ix["gic_log10_kappa_$b"] for b in BLOCKS])]

    ## ---- GlacierMIP3 rungs, GlaMBIE
    rung_y = Dict(b => [Float64(bcrow(b)["com$(replace(string(L), "." => "p"))"]) for L in GMIP_LEVELS] for b in BLOCKS)
    rung_ci = Dict(b => begin
            sig = RUNG_SIG_SCALE .* [Float64(bcrow(b)["sig$(replace(string(L), "." => "p"))"]) for L in GMIP_LEVELS]
            C = (sig * sig') .* (RUNG_CORR .+ (1 - RUNG_CORR) .* Matrix(1.0I, 4, 4))
            inv(C)
        end for b in BLOCKS)
    grate = Dict(b => Float64(bcrow(b).glambie_rate) for b in BLOCKS)

    ## ---- the model: LadrilloModel on 1850-2026, land water :central, tap off; drivers set once; held DAIS parameters
    lm = LadrilloModel(gmst, ohc; years=Y0:Y1, lws=:central, tap=false, ref=FIT_REF)
    m = lm.m
    for (n, s) in HELD
        update_param!(m, _AIS, s, held[n])
    end
    for b in BLOCKS
        update_param!(m, _GLAC, Symbol("glacier_surface_temperature_$b"), tg3[b])
    end
    update_param!(m, _GIS, :greenland_surface_temperature, Float64.(tgis))

    I_ = (smb = [yidx(y) for y in SMB_WIN[1]:SMB_WIN[2]], inv = yidx(INV_YEAR), m1850 = yidx(1850), m1900 = yidx(1900),
          g0 = yidx(GLAMBIE_WIN[1]), g1 = yidx(GLAMBIE_WIN[2]),
          mr0 = yidx(MOUG_REF_WIN[1]), mr1 = yidx(MOUG_REF_WIN[2]), ml0 = yidx(MOUG_LATE_WIN[1]), ml1 = yidx(MOUG_LATE_WIN[2]),
          gisb = [(i0=yidx(w[1]), i1=yidx(w[2]), n=w[2]-w[1]) for w in GISB_WINS])

    return Problem(lm, years, ib, names, [f.comp for f in F], [f.sym for f in F], [f.μ for f in F], [f.σ for f in F],
                   [f.lo for f in F], [f.hi for f in F], np, np + 2*length(SERIES), set_direct, prior_gauss, ix,
                   geo_idx, geo_mu, geo_sd, MvNormal(zeros(length(GEO_SYMS)), Matrix(geo_c)), precip_lo, precip_hi, tbar,
                   gis_tbar, S, B, funch_unit(years), k10_pts,
                   Dict(b => Float64(bcrow(b).S2020_data) for b in BLOCKS), rung_y, rung_ci, grate["R19"],
                   grate["FASTG"] / (grate["SLOWG"] + grate["FASTG"]), held, gap.mean, I_)
end

"""
    logposterior(p, θ) -> Float64

The L27 log-posterior at θ (sampling order, `param_names(p)`): -Inf outside the bounds, the Greenland channel-ordering
wedge or the derived log P0 bounds; otherwise the model is run and every likelihood and prior term summed.
"""
function logposterior(p::Problem, θ::AbstractVector{<:Real})
    np, nk, ix, m = p.np, p.nk, p.ix, p.lm.m
    @inbounds for k in 1:np; (θ[k] < p.lo[k] || θ[k] > p.hi[k]) && return -Inf; end
    (precip_log(p, θ) < p.precip_lo || precip_log(p, θ) > p.precip_hi) && return -Inf
    σn = θ[np+1:2:nk]; ρn = θ[np+2:2:nk]
    (any(σn .<= 0) || any(ρn .< 0) || any(ρn .>= RHO_MAX)) && return -Inf
    ## channel ordering (a wedge in the sampled coordinates): alpha_s <= alpha_f and beta_s <= beta_f
    let r_s = exp(θ[ix["gis_slow_ell"]]), w_s = θ[ix["gis_slow_w"]]
        (w_s * r_s / p.gis_tbar > θ[ix["gis_alpha_f"]]) && return -Inf
        ((1 - w_s) * r_s > θ[ix["gis_beta_f"]])         && return -Inf
    end
    @inbounds for k in p.set_direct
        update_param!(m, p.comp[k], p.sym[k], θ[k])
    end
    update_param!(m, _AIS, :ais_precipitation₀, precip_log(p, θ))
    for b in BLOCKS
        update_param!(m, _GLAC, Symbol("gic_kappa_$b"), 10.0^θ[ix["gic_log10_kappa_$b"]])
    end
    update_param!(m, _GIS, :gis_s_high, 10.0^θ[ix["gis_s_high"]])
    update_param!(m, _AIS, :ais_runoffline_snowheight₀, -θ[ix["ais_runoff_Ton"]] * θ[ix["ais_c"]])
    let r_s = exp(θ[ix["gis_slow_ell"]]), w_s = θ[ix["gis_slow_w"]]
        update_param!(m, _GIS, :gis_alpha_s, w_s * r_s / p.gis_tbar)
        update_param!(m, _GIS, :gis_beta_s, (1 - w_s) * r_s)
    end
    amp = θ[ix["ais_gmst_amp"]]
    update_param!(m, _AIS, :ais_temperature_coefficient, 1.0 / amp)
    update_param!(m, _AIS, :ais_temperature_intercept, -AIS_TANT0 / amp)
    run(m)

    reref(v) = 100 .* (v .- sum(v[p.ib])/length(p.ib))
    Funch = θ[ix["gic_u_unch"]] .* p.funch
    gsic_all_raw  = m[_GLAC, :gsic_sea_level]
    gsic_hind_raw = m[_GLAC, :gsic_hind]
    ais = reref(m[_AIS, :ais_sea_level])
    gsic_flow = reref(gsic_hind_raw .+ Funch)               # hindcast scope: SLOWG + FASTG + F_unch
    gis = reref(m[_GIS, :greenland_sea_level]); te = reref(m[:thermal_expansion, :te_sea_level])
    ll = 0.0
    S = p.series
    for (i, (s, full)) in enumerate(zip((S.ais, S.gsic, S.gis, S.steric), (ais, gsic_flow, gis, te)))
        res = if i == 4
            (full[s.myi] .+ p.d2_steric * [θ[ix["d2_steric_$k"]] for k in 1:D2_BASIS_N]) .- s.obs
        else
            full[s.myi] .- s.obs
        end
        ll += hetero_logl_ar1(res, σn[i], ρn[i], s.ϵ, OBS_CORR_LEN)
    end
    ## Antarctic surface mass balance, 1979-2008 mean (Gt/yr)
    ll += logpdf(Normal(SMB_TARGET_GT, SMB_SIGMA_GT), mean(m[_AIS, :β_total][p.i.smb]) * M3ICE_TO_GT)
    ## glacier inventory: remaining stock in 2000 = sum(a_b) - all-reservoir loss since 1850
    ll += logpdf(Normal(INV_V_M, INV_SIGMA_M),
                 sum(θ[ix["gic_a_$b"]] for b in BLOCKS) - Float64(gsic_all_raw[p.i.inv]))
    ## 1850-1900 loss of SLOWG + FASTG, the two pre-1901 set-asides marginalised over their priors
    let S_1850_1900 = Float64(gsic_hind_raw[p.i.m1900]) - Float64(gsic_hind_raw[p.i.m1850])
        ll += logpdf(Normal(M19_MU_M - LEDGER_SHIFT_M, sqrt(M19_SIGMA_M^2 + LEDGER_SD_M[1]^2 + LEDGER_SD_M[2]^2)),
                     S_1850_1900)
    end
    ## GlacierMIP3 committed loss (%) at four warming levels, per reservoir, in the reservoir's own frame
    for b in BLOCKS
        a = θ[ix["gic_a_$b"]]; bb = θ[ix["gic_b_$b"]]; T0 = θ[ix["gic_T_off_$b"]]
        s20 = p.s2020[b]
        ampb = θ[ix["gic_amp_$b"]]
        r4 = [100.0*(a*(1 - exp(-bb*(ampb*L - T0))) - s20)/max(a - s20, 1e-9) - p.rung_y[b][i]
              for (i, L) in enumerate(GMIP_LEVELS)]
        ll += -0.5 * (r4' * (p.rung_ci[b] * r4))
    end
    ## GlaMBIE 2000-2024: the FASTG share of SLOWG + FASTG
    let n = GLAMBIE_WIN[2] - GLAMBIE_WIN[1],
        rs = 1000.0*(Float64(m[_GLAC, :gsic_slowg][p.i.g1]) - Float64(m[_GLAC, :gsic_slowg][p.i.g0])) / n,
        rf = 1000.0*(Float64(m[_GLAC, :gsic_fastg][p.i.g1]) - Float64(m[_GLAC, :gsic_fastg][p.i.g0])) / n
        tot = rs + rf
        abs(tot) > GLAMBIE_TOT_FLOOR && (ll += logpdf(Normal(p.glambie_fast_share, GLAMBIE_SHARE_SD), rf / tot))
    end
    ## GlaMBIE 2000-2024 R19 rate (mm/yr)
    let r19rate = 1000.0*(Float64(m[_GLAC, :gsic_r19][p.i.g1]) - Float64(m[_GLAC, :gsic_r19][p.i.g0])) /
                  (GLAMBIE_WIN[2] - GLAMBIE_WIN[1])
        ll += logpdf(Normal(p.r19_rate_mu, R19_RATE_SD), r19rate)
    end
    ## Mouginot 2019: surface (fast) share of the extra loss rate of 2000-2018 over 1972-1990
    let tot = m[_GIS, :greenland_sea_level], fst = m[_GIS, :gis_fast], q = p.i
        rate(x, i0, i1, n) = (Float64(x[i1]) - Float64(x[i0])) / n
        nref, nlate = MOUG_REF_WIN[2]-MOUG_REF_WIN[1], MOUG_LATE_WIN[2]-MOUG_LATE_WIN[1]
        d_tot  = rate(tot, q.ml0, q.ml1, nlate) - rate(tot, q.mr0, q.mr1, nref)
        d_fast = rate(fst, q.ml0, q.ml1, nlate) - rate(fst, q.mr0, q.mr1, nref)
        abs(d_tot) > 1e-12 && (ll += logpdf(Normal(MOUG_SHARE, MOUG_SHARE_SD), d_fast / d_tot))
    end
    ## Mouginot 2019: the active basin's share of the mean loss rate in each window (high follows by sum-to-one)
    let fa = m[_GIS, :gis_fast_active], sa = m[_GIS, :gis_slow_active], fh = m[_GIS, :gis_fast_high],
        sh = m[_GIS, :gis_slow_high]
        sl(f, s, i) = Float64(f[i]) + Float64(s[i])
        for (wi, tgt) in zip(p.i.gisb, GISB_SHARE_ACTIVE)
            d_a = (sl(fa, sa, wi.i1) - sl(fa, sa, wi.i0)) / wi.n
            d_h = (sl(fh, sh, wi.i1) - sl(fh, sh, wi.i0)) / wi.n
            dtot = d_a + d_h
            abs(dtot) > GISB_TOT_FLOOR && (ll += logpdf(Normal(tgt, GISB_SHARE_SD), d_a / dtot))
        end
    end
    ## priors
    lp = 0.0
    @inbounds for k in p.prior_gauss
        lp += logpdf(Normal(p.mu[k], p.sigma[k]), θ[k])
    end
    for b in BLOCKS
        lp += logpdf(Normal(k10c(p, b, θ[ix["gic_amp_$b"]]), K10_SIG), θ[ix["gic_log10_kappa_$b"]])
    end
    let g = θ[p.geo_idx]
        g[findfirst(==(ix["ais_precip_u"]), p.geo_idx)] = precip_log(p, θ)
        lp += logpdf(p.geo_prior, (g .- p.geo_mu) ./ p.geo_sd)
    end
    for i in 1:length(SERIES); lp += logpdf(truncated(Normal(0, NOISE_SD_PRIOR), 0, Inf), σn[i]); end
    return ll + lp
end

"""
    start_point(p, seed; file) -> Vector{Float64}

The chain's start: row i of the over-dispersed starts file for seed SEEDS_CANON[i] (or SEEDS_REPLICATE[i], the same
row on another RNG stream). The rows are posterior draws of a tuning run at ais_iceflow0 quantiles 0.02/0.35/0.65/0.98;
columns are read by name, so the file may carry parameters this calibration does not sample.
"""
function start_point(p::Problem, seed::Integer; file=datapath("calibration", STARTS_FILE))
    si = something(findfirst(==(seed), SEEDS_CANON), findfirst(==(seed), SEEDS_REPLICATE), 0)
    si == 0 && error("no start row for seed $seed (canonical $SEEDS_CANON, replicate $SEEDS_REPLICATE)")
    st = CSV.read(file, DataFrame)
    nrow(st) >= si || error("$file has $(nrow(st)) rows, need >= $si")
    missing_cols = [n for n in p.names if !hasproperty(st, Symbol(n))]
    isempty(missing_cols) || error("$file lacks $(join(missing_cols, ", "))")
    return [Float64(st[si, Symbol(n)]) for n in p.names]
end

"""
    proposal_covariance(p; file) -> (cov0, message)

The sampler's initial proposal covariance: a diagonal default (0.1 x min(prior sd, range/4); 0.02 x the paleo sd for the
geometry block; 0.3 and 0.1 for the noise sd and rho), overwritten BY NAME with the rows and columns of an adapted
covariance from an earlier run (data/calibration/adapted_cov_L26_named.csv). Checks positive definiteness and the
GEO_SEED_FLOOR sd floors.
"""
function proposal_covariance(p::Problem; file=datapath("calibration", ADCOV_FILE))
    prop = vcat([PROP_SCALE*Float64(min(p.sigma[k], (p.hi[k] - p.lo[k])/4)) for k in 1:p.np],
                repeat(collect(NOISE_PROP), length(SERIES)))
    for k in p.geo_idx; prop[k] = GEO_PROP_SCALE * Float64(p.sigma[k]); end
    cov0 = Matrix(Diagonal(prop.^2))
    adf = CSV.read(file, DataFrame)
    size(adf, 1) == size(adf, 2) || error("$file is not square")
    old, onames = Matrix(adf), names(adf)
    oi, ni = Int[], Int[]
    for (i, nm) in enumerate(onames)
        j = findfirst(==(nm), p.names)
        isnothing(j) && continue
        push!(oi, i); push!(ni, j)
    end
    cov0[ni, ni] = old[oi, oi]
    dropped, fresh = setdiff(onames, p.names), setdiff(p.names, onames)
    msg = "proposal seed: $(length(oi)) of $(size(old, 1)) rows of $(basename(file)) mapped by name" *
          (isempty(dropped) ? "" : "; not sampled here: " * join(dropped, ", ")) *
          (isempty(fresh) ? "" : "; diagonal default for " * join(fresh, ", "))
    isposdef(cov0) || error("seed proposal covariance is not positive definite")
    sd = sqrt.(diag(cov0))
    gated = vcat([n == "ais_precip0_LOG" ? "ais_precip_u" : n for n in GEO_NAMES], ["gis_s_high"])
    bad = [(n, sd[p.ix[n]]) for n in gated if sd[p.ix[n]] < GEO_SEED_FLOOR[n == "ais_precip_u" ? "ais_precip0_LOG" : n]]
    isempty(bad) || error("proposal seed sd below its floor (the parameter would be frozen): $bad")
    return cov0, msg
end

"""
    parameter_table(p) -> DataFrame

One row per sampled parameter in sampling order, with the prior as `logposterior` scores it (cf. the research
calibrator's --dump-priors, data/constants/ladrillo_priors_L27.csv, which is Table A1's input).
"""
function parameter_table(p::Problem)
    rows = DataFrame(index=Int[], name=String[], comp=String[], mu=Float64[], sigma=Float64[], lo=Float64[],
                     hi=Float64[], prior_form=String[])
    kap = Dict(p.ix["gic_log10_kappa_$b"] => b for b in BLOCKS)
    for k in 1:p.np
        form = k in p.geo_idx ? "joint paleo MvNormal" :
               haskey(kap, k) ? "N(k10c(gic_amp_$(kap[k])), $K10_SIG)" :
               p.sigma[k] >= 10.0 ? "flat" : "N(mu, sigma)"
        push!(rows, (k, p.names[k], String(p.comp[k]), p.mu[k], p.sigma[k], p.lo[k], p.hi[k], form))
    end
    for (i, s) in enumerate(SERIES)
        push!(rows, (p.np + 2i - 1, "sd_$s", "noise", 0.0, NOISE_SD_PRIOR, 0.0, Inf, "half-normal"))
        push!(rows, (p.np + 2i, "rho_$s", "noise", NaN, NaN, 0.0, RHO_MAX, "flat"))
    end
    return rows
end

"""
    run_chain(p, n, seed; show_progress=true) -> NamedTuple

One chain of `n` iterations of the robust adaptive Metropolis sampler (Vihola 2012; RobustAdaptiveMetropolisSampler.jl,
target acceptance OPT_ALPHA) from `start_point(p, seed)` with the `proposal_covariance(p)` seed. The global RNG is seeded
with `seed` immediately before sampling and nothing else draws from it. Adaptation does not depend on `n`, so a shorter
run is a prefix of a longer one with the same seed. Returns the chain (n x nk), its log-posterior, the acceptance rate,
the adapted covariance, the start and its log-posterior, and the proposal-seed message.
"""
function run_chain(p::Problem, n::Integer, seed::Integer; show_progress::Bool=true)
    θ0 = start_point(p, seed)
    lp0 = logposterior(p, θ0)
    isfinite(lp0) || error("logposterior at the start row for seed $seed is $lp0: every Metropolis ratio would be NaN")
    cov0, msg = proposal_covariance(p)
    Random.seed!(seed)
    r = RAM_sample(θ -> logposterior(p, θ), θ0, cov0, Int(n); opt_α=OPT_ALPHA, output_log_probability_x=true,
                   show_progress=show_progress)
    return (chain=r.chain, logpost=r.log_probabilities_x, accept=r.acceptance_rate, cov=r.M, θ0=θ0, lp0=lp0,
            seed_message=msg)
end

end # module Calibration
