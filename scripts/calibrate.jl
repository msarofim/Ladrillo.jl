# Sample the L27 calibration (src/calibration.jl): one chain per call.
#
#   julia --project=. --threads=1 scripts/calibrate.jl <n_iter> <seed> [--outdir=outputs/calibration]
#
# Production: n_iter = 2,000,000 and seeds 2026, 2027, 2028, 2029, one process per seed (about 4 h each on an Apple M4,
# ~2 GB of memory and ~2 GB of output each). Seed i starts from row i of data/calibration/overdispersed_starts_L27.csv;
# seeds 3026-3029 run the same start rows on other RNG streams. scripts/postprocess_calibration.jl (not yet ported)
# pools the four chains into the posterior.
#
# Writes to <outdir>:
#   chain_L27_seed<seed>_n<n_iter>.csv              the parameters in sampling order, log_post, accept_rate (one row per
#                                                   iteration)
#   adapted_cov_L27_seed<seed>_n<n_iter>.csv        the sampler's final proposal covariance, header = parameter names
#   chain_L27_seed<seed>_n<n_iter>_provenance.txt   how the chain was made, including the seed and where it is applied
#
# REPRODUCIBILITY. The global RNG is seeded with <seed> immediately before sampling and nothing else draws from it; BLAS
# runs on one thread (set below). The sampler's adaptation does not depend on n_iter, so a shorter run is a prefix of a
# longer one. On the reference runtime (Apple M4, Julia 1.12.6, the Manifest's versions) the chain is bit-identical to
# the research production chain (test/calibration/gate_chain_prefix.jl); elsewhere expect statistical agreement only.

using Ladrillo, CSV, DataFrames, Printf
import LinearAlgebra

const C = Ladrillo.Calibration
arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const POS = filter(a -> !startswith(a, "--"), ARGS)
length(POS) == 2 || error("usage: scripts/calibrate.jl <n_iter> <seed> [--outdir=...]")
const N_ITER, SEED = parse(Int, POS[1]), parse(Int, POS[2])
const OUTDIR = arg("--outdir=", joinpath(dirname(@__DIR__), "outputs", "calibration"))
const BLAS_THREADS = 1
LinearAlgebra.BLAS.set_num_threads(BLAS_THREADS)
mkpath(OUTDIR)
const STEM = "chain_$(C.TAG)_seed$(SEED)_n$(N_ITER)"

p = C.problem()
println("calibration $(C.TAG): $(p.np) physical + $(p.nk - p.np) noise parameters; held (not sampled): " *
        join(["$n = $(p.held[n])" for (n, _) in C.HELD], ", "))
t = @elapsed r = C.run_chain(p, N_ITER, SEED)
@printf("logposterior(start) = %s | %d iterations in %.0f s | acceptance %.4f\n", repr(r.lp0), N_ITER, t, r.accept)
println(r.seed_message)

df = DataFrame(r.chain, C.param_names(p))
df.log_post = r.logpost
df.accept_rate = fill(r.accept, nrow(df))
CSV.write(joinpath(OUTDIR, STEM * ".csv"), df)
CSV.write(joinpath(OUTDIR, "adapted_cov_$(C.TAG)_seed$(SEED)_n$(N_ITER).csv"), DataFrame(r.cov, C.param_names(p)))
si = something(findfirst(==(SEED), C.SEEDS_CANON), findfirst(==(SEED), C.SEEDS_REPLICATE))
open(joinpath(OUTDIR, STEM * "_provenance.txt"), "w") do io
    println(io, "Ladrillo $(pkgversion(Ladrillo)) | scripts/calibrate.jl $N_ITER $SEED | calibration $(C.TAG)")
    println(io, "seed $SEED: Random.seed!($SEED) on the global RNG immediately before RAM_sample (nothing else draws " *
                "from it); start = row $si of data/calibration/$(C.STARTS_FILE); logposterior(start) = $(repr(r.lp0))")
    println(io, "sampler: RobustAdaptiveMetropolisSampler $(pkgversion(parentmodule(C.RAM_sample))) RAM_sample, " *
                "opt_alpha $(C.OPT_ALPHA), $N_ITER iterations, acceptance $(r.accept)")
    println(io, r.seed_message)
    println(io, "held (not sampled): " * join(["$n = $(p.held[n])" for (n, _) in C.HELD], ", ") *
                " (data/constants/central_row.csv; see Calibration.HELD)")
    println(io, "forcing: data/forcing/fair_mean_{gmst,ohc}_$(C.FORCING).csv, FaIR 2.2.4 (calib 1.6.0) ensemble mean " *
                "($(Ladrillo.fair_seed_note(C.FORCING))) | model $(C.Y0)-$(C.Y1), land water :central, Greenland tap off")
    println(io, "runtime: julia $(VERSION), $(Sys.MACHINE), BLAS threads $(LinearAlgebra.BLAS.get_num_threads()), " *
                "julia threads $(Threads.nthreads()), $(round(t, digits=1)) s")
end
println("wrote $(joinpath(OUTDIR, STEM)).csv (+ adapted covariance, provenance)")
