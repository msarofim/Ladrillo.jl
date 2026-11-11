# Post-process the L27 calibration chains into the posterior file (src/postprocess.jl).
#
#   julia --project=. --threads=1 scripts/postprocess_calibration.jl [--chains=DIR] [--outdir=DIR] [--accept-slr] [--force]
#
# Production (the run behind data/posterior/ladrillo_L27_posterior.csv): the four 2,000,000-iteration chains of seeds
# 2026-2029 from scripts/calibrate.jl, with --accept-slr. About 25 min on an Apple M4 (21 of them the parameter ESS),
# 3.3 GB peak memory.
#
# --chains   directory holding chain_L27_seed<seed>_n<n_iter>.csv for each seed in 2026-2029 (default outputs/calibration)
# --outdir   where to write (default outputs/calibration)
# --accept-slr  write the posterior when the projected sea-level certificate passes even if some parameter marginals
#               do not (the criterion L27 was accepted on)
# --force    write anyway, with a _NOTCONVERGED suffix
#
# Writes to <outdir>:
#   slr_convergence_L27.csv          the certificate: R-hat and ESS of projected total sea level at 2100 and 2150
#   parameter_convergence_L27.csv    R-hat, ESS and tau of every parameter and log_post
#   ladrillo_L27_posterior.csv       the 10,000-draw posterior (the format of data/posterior/), if accepted
#   adapted_cov_L27.csv              covariance of the pooled draws, header = parameter names, if accepted
#   postprocess_L27_provenance.txt   inputs, settings, the seed and where it is applied, the criterion that passed
# Nothing here draws random numbers except the certificate's paleo attachment (fixed seed, see the provenance file).

using Ladrillo, CSV, DataFrames, Printf, Statistics
import LinearAlgebra

const P = Ladrillo.Postprocess
arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const ROOT = dirname(@__DIR__)
const CHAINDIR = arg("--chains=", joinpath(ROOT, "outputs", "calibration"))
const OUTDIR = arg("--outdir=", joinpath(ROOT, "outputs", "calibration"))
const ACCEPT_SLR = "--accept-slr" in ARGS
const FORCE = "--force" in ARGS
const BLAS_THREADS = 1
LinearAlgebra.BLAS.set_num_threads(BLAS_THREADS)
mkpath(OUTDIR)
t_start = time()

## 1. read
files, n_iter = P.chain_files(CHAINDIR)
labels = ["seed$s" for s in P.SEEDS]
println("Post-processing calibration $(P.TAG): $(length(files)) chains of $n_iter iterations from $CHAINDIR")
chains = P.read_chains(files)
const N_POOL = sum(nrow.(chains))
const STEP = P.subsample_step(N_POOL)
for (f, c) in zip(files, chains)
    @printf("  %s  (%d post-burn)\n", basename(f), nrow(c))
end

## 2. noise mode
println("\n=== NOISE-MODE CHECK: second-half median sd_gis < $(P.NOISE_SD_GIS_MAX) cm in every chain ===")
nm = P.noise_mode(chains)
for (lab, r) in zip(labels, eachrow(nm))
    @printf("  %s: sd_gis median %.4f  sd_ais %.4f  log_post %.1f  %s\n", lab, r.sd_gis, r.sd_ais, r.log_post,
            r.pass ? "OK" : "*** NOISE MODE ***")
end
all(nm.pass) || error("noise-mode check failed on $(labels[.!nm.pass]): not post-processing")
println("NOISE-MODE CHECK: PASS")

## 3. certificate (printed before the parameters, in the research run's order)
@printf("\nProjected sea-level certificate: %s FaIR mean, %d-%d, land water :%s, tap %s, %d draws per chain\n",
        P.CERT_SCENARIO, first(P.CERT_YEARS), last(P.CERT_YEARS), P.CERT_LWS, P.CERT_TAP ? "on" : "off", P.CERT_PER_CHAIN)
t0 = time()
cert, totals = P.slr_certificate(chains, n_iter)
@printf("  %d model runs in %.0f s\n", length(chains) * P.CERT_PER_CHAIN, time() - t0)
q(v, p) = quantile(v, p)
for y in P.CERT_HORIZONS
    @printf("\nSLR @%d (cm, rel. %d-%d)\n", y, Ladrillo.REF_PERIOD[1], Ladrillo.REF_PERIOD[2])
    @printf("  %-10s %8s %8s %8s %8s %8s\n", "chain", "q05", "q50", "q95", "mean", "sd")
    for (ci, lab) in enumerate(labels)
        v = totals[y][:, ci]
        @printf("  %-10s %8.2f %8.2f %8.2f %8.2f %8.2f\n", lab, q(v, 0.05), q(v, 0.50), q(v, 0.95), mean(v), std(v))
    end
    v = vec(totals[y])
    @printf("  %-10s %8.2f %8.2f %8.2f %8.2f %8.2f\n", "POOLED", q(v, 0.05), q(v, 0.50), q(v, 0.95), mean(v), std(v))
end
println("\nCROSS-CHAIN CONVERGENCE OF PROJECTED SLR ($(length(chains)) chains x $(P.CERT_PER_CHAIN) draws)")
@printf("  %-14s %8s %10s %12s %12s %10s\n", "quantity", "R-hat", "ESS", "sd(medians)", "mean(sd_wc)", "ratio")
for r in eachrow(cert)
    @printf("  SLR@%-10d %8.3f %10.1f %12.3f %12.3f %10.3f\n", r.horizon, r.rhat, r.ess, r.sd_medians, r.mean_sd_wc,
            r.sd_medians / r.mean_sd_wc)
end
CSV.write(joinpath(OUTDIR, "slr_convergence_$(P.TAG).csv"), cert)
const SLR_OK = all(isfinite, cert.rhat) && all(cert.rhat .< P.RHAT_MAX)
println("certificate: projected SLR is $(SLR_OK ? "" : "NOT ")converged across chains (test: R-hat < $(P.RHAT_MAX) at every horizon)")

## 4. parameters
t0 = time()
pc = P.parameter_convergence(chains)
println("\n$(length(chains)) chains × $(nrow(chains[1])) draws. Convergence (target R̂<$(P.RHAT_MAX), ESS>$(P.ESS_MIN)):")
for r in eachrow(pc)
    r.pass || @printf("  %-24s R̂=%.3f ESS=%.1f τ=%.0f  <-- check\n", r.name, r.rhat, r.ess, r.tau)
end
const CONVERGED = all(pc.pass)
CONVERGED ? println("  all params converged (R̂<$(P.RHAT_MAX), ESS>$(P.ESS_MIN)).") :
            println("  $(count(.!pc.pass)) params NOT converged.")
@printf("  (%.0f s)\n", time() - t0)
CSV.write(joinpath(OUTDIR, "parameter_convergence_$(P.TAG).csv"), pc)

## 5. acceptance
slr_accepted = !CONVERGED && ACCEPT_SLR && SLR_OK
if !CONVERGED && ACCEPT_SLR
    println("\n--accept-slr: deliverable-level convergence:")
    for r in eachrow(cert)
        @printf("  SLR@%d  R̂=%.3f  ESS=%.1f\n", r.horizon, r.rhat, r.ess)
    end
    println(slr_accepted ? "ACCEPTED ON THE CERTIFICATE: parameter marginals not converged, projected SLR R̂<$(P.RHAT_MAX) " *
                           "at all horizons." : "certificate R̂ >= $(P.RHAT_MAX) at some horizon: NOT accepted.")
end
criterion = CONVERGED ? "every parameter (R-hat <= $(P.RHAT_MAX), ESS >= $(P.ESS_MIN))" :
            slr_accepted ? "the projected sea-level certificate (--accept-slr); $(count(.!pc.pass)) parameter marginals " *
                           "not converged: $(join(pc.name[.!pc.pass], ", "))" :
            FORCE ? "none (--force): written with the _NOTCONVERGED suffix" : "none"

## 6-7. subsample and covariance
written = String[]
if CONVERGED || slr_accepted || FORCE
    sfx = (CONVERGED || slr_accepted) ? "" : "_NOTCONVERGED"
    draws, C = P.subsample(chains)
    chains = nothing; GC.gc()
    out = joinpath(OUTDIR, "ladrillo_$(P.TAG)_posterior$(sfx).csv")
    CSV.write(out, draws)
    @printf("\nWrote %s  (%d draws: every %d-th of %d pooled)\n", out, nrow(draws), STEP, N_POOL)
    cout = joinpath(OUTDIR, "adapted_cov_$(P.TAG)$(sfx).csv")
    CSV.write(cout, DataFrame(C, names(draws)))
    println("Wrote $cout (covariance of the pooled draws)")
    push!(written, basename(out), basename(cout))
else
    println("\nNOT CONVERGED: no posterior written. Rerun longer, pass --accept-slr to accept on the certificate, or " *
            "--force to write with a _NOTCONVERGED suffix.")
end

## provenance
open(joinpath(OUTDIR, "postprocess_$(P.TAG)_provenance.txt"), "w") do io
    println(io, "Ladrillo $(pkgversion(Ladrillo)) | scripts/postprocess_calibration.jl $(join(ARGS, " ")) | calibration $(P.TAG)")
    println(io, "chains, pooled in this order (second half of each):")
    for f in files
        println(io, "  $f ($(filesize(f)) bytes, $n_iter iterations)")
    end
    println(io, "noise-mode check: second-half median sd_gis " * join([@sprintf("%.4f", x) for x in nm.sd_gis], ", ") *
                " cm, all < $(P.NOISE_SD_GIS_MAX)")
    println(io, "certificate (slr_convergence_$(P.TAG).csv): total sea level at $(join(P.CERT_HORIZONS, ", ")), rel. " *
                "$(Ladrillo.REF_PERIOD[1])-$(Ladrillo.REF_PERIOD[2]); $(P.CERT_PER_CHAIN) draws per chain evenly spaced " *
                "through its second half; $(P.CERT_SCENARIO) FaIR 2.2.4 (calib 1.6.0) ensemble mean " *
                "($(Ladrillo.fair_seed_note(P.CERT_SCENARIO))), $(first(P.CERT_YEARS))-$(last(P.CERT_YEARS)), " *
                "land water :$(P.CERT_LWS), Greenland tap $(P.CERT_TAP ? "on" : "off")")
    println(io, "seed $(P.PALEO_SEED): MersenneTwister($(P.PALEO_SEED)) in attach_paleo!, once per chain, picks the " *
                "$(P.CERT_PER_CHAIN) (lambda, T_crit) rows of data/constants/paleo_fastdyn_draws.csv given to that " *
                "chain's certificate draws (the same rows for every chain); gamma = its paleo median. No other random numbers.")
    println(io, "parameters (parameter_convergence_$(P.TAG).csv): MCMCDiagnosticTools $(pkgversion(parentmodule(P.rhat))) " *
                "rhat and ess, ESS maxlag min(draws - 4, $(P.ESS_MAXLAG)); pass = R-hat <= $(P.RHAT_MAX) and ESS >= $(P.ESS_MIN)")
    println(io, "accepted on: $criterion")
    isempty(written) || println(io, "posterior: every $(STEP)-th of the $(N_POOL) pooled draws, at most " *
                                    "$(P.N_SUB); covariance of all pooled draws + $(P.COV_JITTER) I. Wrote " *
                                    join(written, ", "))
    println(io, "runtime: julia $(VERSION), $(Sys.MACHINE), BLAS threads $(LinearAlgebra.BLAS.get_num_threads()), " *
                "julia threads $(Threads.nthreads()), $(round(time() - t_start, digits=1)) s")
end
println("wrote $(joinpath(OUTDIR, "postprocess_$(P.TAG)_provenance.txt"))")
