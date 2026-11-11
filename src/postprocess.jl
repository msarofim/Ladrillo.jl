# Post-processing of the calibration (L27): from the four chains to the posterior file.
#
# A port of the research post-processing of 2026-09-20 (SLR-RFF-BRICK run_l27_postprocess.sh): its noise-mode check,
# julia/diag_slr_convergence_by_chain_ladrillo.jl (the projected sea-level convergence certificate) and
# julia/postprocess_mcmc_ext.jl --tag=L27 --accept-slr (parameter convergence, pooling, subsample, pooled covariance).
# scripts/postprocess_calibration.jl runs it; test/calibration/gate_postprocess.jl checks it against the research run.
#
# STEPS
#   1. read          each chain file of the canonical seeds, in seed order; keep its second half.
#   2. noise mode    the median sd_gis of each second half must be below NOISE_SD_GIS_MAX. A chain caught in the noise
#                    mode (the Greenland AR(1) sd absorbing the misfit) sat at 0.6-1.0 cm in the L27 tuning run.
#   3. parameters    R-hat and ESS of every parameter and of log_post across the chains.
#   4. certificate   R-hat and ESS of projected total sea level at CERT_HORIZONS across the chains, from CERT_PER_CHAIN
#                    draws evenly spaced through each second half.
#   5. acceptance    the posterior is written when every parameter passes (R-hat <= RHAT_MAX and ESS >= ESS_MIN), or,
#                    with `accept_slr`, when the certificate passes at every horizon. The AIS geometry parameters lie on
#                    a ridge whose marginals do not converge in feasible compute while the projection does (criterion
#                    adopted 2026-07-19). L27 was accepted on the certificate: 8 marginals did not pass.
#   6. subsample     pool the second halves in seed order and take every (pooled / N_SUB)-th draw: N_SUB draws, rows
#                    1-2500 from seed 2026. No random numbers.
#   7. covariance    of the pooled draws plus COV_JITTER * I: a proposal seed for a later calibration.
#
# THE CERTIFICATE'S MODEL is the projection model as it stood when L27 was accepted: SSP2-4.5 FaIR ensemble-mean
# forcing, 1850-2150, land water :central, Greenland tap off, sea level relative to REF_PERIOD, (lambda, T_crit) attached
# as joint paleo draws in one block per chain (seed PALEO_SEED; the four chains share the same rows) and gamma at its
# paleo median. The shipped projections differ in land water (:observed) and in the paleo assignment; neither enters
# a cross-chain comparison of one chain against another.

module Postprocess

using CSV, DataFrames, Statistics
import LinearAlgebra: I
import MCMCDiagnosticTools: rhat, ess
using ..Ladrillo: LadrilloModel, mean_forcing, run_draw!, sea_level, attach_paleo!, native_greenland!, REF_PERIOD,
                  PALEO_SEED
import ..Calibration

const TAG = Calibration.TAG
const SEEDS = Calibration.SEEDS_CANON
const RHAT_MAX = 1.05            # parameters and certificate alike
const ESS_MIN = 400
## ESS needs an explicit maxlag. MCMCDiagnosticTools' default (250) floors the estimate at n/500, which over-reported the
## AIS geometry block's ESS 150-270 fold; maxlag = n returns NaN on chains of 1e6 draws or more. Capped here, and a
## non-finite ESS fails.
const ESS_MAXLAG = 200_000
const N_SUB = 10_000
const COV_JITTER = 1e-10
const NOISE_SD_GIS_MAX = 0.10    # cm
const CERT_SCENARIO = "ssp245"
const CERT_YEARS = 1850:2150
const CERT_HORIZONS = (2100, 2150)
const CERT_PER_CHAIN = 400
const CERT_LWS = :central
const CERT_TAP = false

"""The research chains name two glacier reservoirs SLOWP and FAST; the package names them SLOWG and FASTG."""
package_name(s::AbstractString) = replace(replace(s, r"_SLOWP$" => "_SLOWG"), r"_FAST$" => "_FASTG")

const CHAIN_RX = Regex("^chain_" * TAG * "_seed(\\d+)_n(\\d+)\\.csv\$")

"""
    chain_files(dir; seeds=SEEDS) -> (files, n_iter)

The chain file of each seed in `dir` (`chain_L27_seed<seed>_n<n_iter>.csv`, as scripts/calibrate.jl writes them), in
seed order, which is the order they are pooled in. Every seed needs exactly one file and all files the same n_iter.
"""
function chain_files(dir::AbstractString; seeds=SEEDS)
    found = [(f, parse(Int, m[1]), parse(Int, m[2])) for (f, m) in ((f, match(CHAIN_RX, f)) for f in readdir(dir))
             if m !== nothing]
    files, ns = String[], Int[]
    for s in seeds
        hit = [(f, n) for (f, sd, n) in found if sd == s]
        length(hit) == 1 || error("seed $s: $(length(hit)) chain files in $dir ($(join(first.(hit), ", ")))")
        push!(files, joinpath(dir, hit[1][1])); push!(ns, hit[1][2])
    end
    allequal(ns) || error("chains of different lengths: " * join(["$(basename(f)) (n $n)" for (f, n) in zip(files, ns)], ", "))
    other = [f for (f, sd, _) in found if !(sd in seeds)]
    isempty(other) || @info "not pooled (seed not in $(collect(seeds))): $(join(other, ", "))"
    return files, ns[1]
end

"""
    read_chains(files) -> Vector{DataFrame}

The second half of each chain: the 50 parameters in sampling order (package names) and log_post. The header must be
`Calibration.param_names` (research or package glacier names), log_post, accept_rate.
"""
function read_chains(files)
    expected = vcat(Calibration.param_names(Calibration.problem()), ["log_post", "accept_rate"])
    chains = DataFrame[]
    for f in files
        hdr = String.(propertynames(CSV.read(f, DataFrame; limit=0)))
        package_name.(hdr) == expected || error("$(basename(f)): header is not the L27 parameter list " *
                                                "(differs at $(join(symdiff(package_name.(hdr), expected), ", ")))")
        d = CSV.read(f, DataFrame; select=hdr[1:end-1])
        rename!(d, package_name.(names(d)))
        push!(chains, d[(nrow(d) ÷ 2 + 1):end, :])
        d = nothing; GC.gc()
    end
    allequal(nrow.(chains)) || error("second halves of different lengths: $(nrow.(chains))")
    return chains
end

"""Median sd_gis, sd_ais and log_post of each chain's second half, and whether sd_gis is below NOISE_SD_GIS_MAX."""
noise_mode(chains) = DataFrame(chain=eachindex(chains), sd_gis=[median(c.sd_gis) for c in chains],
                               sd_ais=[median(c.sd_ais) for c in chains], log_post=[median(c.log_post) for c in chains],
                               pass=[median(c.sd_gis) < NOISE_SD_GIS_MAX for c in chains])

"""
    parameter_convergence(chains) -> DataFrame

R-hat and ESS across the chains of every parameter and of log_post (MCMCDiagnosticTools defaults, ESS maxlag capped at
ESS_MAXLAG), the implied autocorrelation time tau = draws / ESS, and whether R-hat <= RHAT_MAX and ESS >= ESS_MIN.
"""
function parameter_convergence(chains)
    n, nc = nrow(chains[1]), length(chains)
    out = DataFrame(name=String[], rhat=Float64[], ess=Float64[], tau=Float64[], pass=Bool[])
    for p in names(chains[1])
        arr = Array{Float64}(undef, n, nc)
        for (ci, ch) in enumerate(chains); arr[:, ci] = Float64.(ch[1:n, p]); end
        r = rhat(arr); e = ess(arr; maxlag=min(n - 4, ESS_MAXLAG))
        push!(out, (p, r, e, (n * nc) / max(e, 1e-9),
                    isfinite(r) && isfinite(e) && r <= RHAT_MAX && e >= ESS_MIN))
    end
    return out
end

"""
    slr_certificate(chains, n_iter) -> (table, totals)

Projected total sea level (cm rel. REF_PERIOD) at CERT_HORIZONS for CERT_PER_CHAIN draws per chain, evenly spaced through
each second half, on the certificate's model (header comment). `table` has one row per horizon: R-hat and ESS across
the chains (ESS maxlag = draws per chain), the sd of the chain medians, the mean within-chain sd, and the run's shape.
`totals[y]` is a draws x chains matrix.
"""
function slr_certificate(chains, n_iter::Integer)
    years = collect(CERT_YEARS)
    g, o = mean_forcing(CERT_SCENARIO; years=years)
    lm = LadrilloModel(g, o; years=years, lws=CERT_LWS, tap=CERT_TAP, ref=REF_PERIOD)
    iy = [findfirst(==(y), years) for y in CERT_HORIZONS]
    totals = Dict(y => Matrix{Float64}(undef, CERT_PER_CHAIN, length(chains)) for y in CERT_HORIZONS)
    for (ci, ch) in enumerate(chains)
        step = max(1, nrow(ch) ÷ CERT_PER_CHAIN)
        rows = collect(1:step:nrow(ch))
        length(rows) >= CERT_PER_CHAIN || error("chain $ci: $(nrow(ch)) draws give $(length(rows)) at step $step")
        d = ch[rows[1:CERT_PER_CHAIN], :]
        attach_paleo!(d)                  # one block over this chain's draws
        native_greenland!(d)
        for (j, r) in enumerate(eachrow(d))
            run_draw!(lm, r)
            tot = sea_level(lm, :total)
            for (k, y) in enumerate(CERT_HORIZONS); totals[y][j, ci] = tot[iy[k]]; end
        end
    end
    table = DataFrame(horizon=Int[], rhat=Float64[], ess=Float64[], sd_medians=Float64[], mean_sd_wc=Float64[],
                      n_chains=Int[], n_per_chain=Int[], niter=Int[], gis_shape=Bool[])
    for y in CERT_HORIZONS
        arr = totals[y]
        cols = [arr[:, ci] for ci in axes(arr, 2)]
        push!(table, (y, rhat(arr), ess(arr; maxlag=size(arr, 1)), std(median.(cols)), mean(std.(cols)),
                      size(arr, 2), size(arr, 1), n_iter, true))   # gis_shape: the amp(GMST) law is always on here
    end
    return table, totals
end

"""Thinning step of the subsample for `n` pooled draws."""
subsample_step(n::Integer) = max(1, n ÷ N_SUB)

"""
    subsample(chains) -> (draws, pooled_cov)

Pool the second halves in chain order, take every (pooled / N_SUB)-th row (N_SUB rows, no random numbers), drop
log_post. `pooled_cov` is the covariance of all pooled draws plus COV_JITTER * I.
"""
function subsample(chains)
    pool = vcat(chains...)
    step = subsample_step(nrow(pool))
    pars = [p for p in names(pool) if p != "log_post"]
    draws = pool[1:step:end, :][1:min(N_SUB, end), pars]
    M = Matrix{Float64}(pool[:, pars])
    pool = nothing; GC.gc()
    return draws, cov(M) .+ COV_JITTER * I(size(M, 2))
end

end # module Postprocess
