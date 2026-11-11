# Post-processing gate: does scripts/postprocess_calibration.jl turn the production L27 chains into the shipped posterior?
#
#   julia --project=. --threads=1 test/calibration/gate_postprocess.jl <SLR-RFF-BRICK repo> [--outdir=DIR] [--reuse]
#
# Exit 0 = PASS. Checks the chains' sha256 against test/reference/SOURCES.csv, runs the script on them with --accept-slr
# (about 25 min, 3.3 GB peak memory; output and log in --outdir, default outputs/calibration/gate_postprocess), then requires
#   (1) ladrillo_L27_posterior.csv byte-identical to data/posterior/ladrillo_L27_posterior.csv (the research subsample of
#       2026-09-20 with SLOWP/FAST renamed SLOWG/FASTG)
#   (2) slr_convergence_L27.csv byte-identical to the research certificate of 2026-09-20
#   (3) adapted_cov_L27.csv: header = the posterior's header, every value == the research pooled covariance (whose header
#       is x1..x50)
#   (4) the run's log lines matching LOG_PATTERNS == the research log's (test/reference/research_postprocess_L27_log_lines.txt):
#       the noise check, the certificate's per-chain quantiles and R-hat, the 8 parameters that fail with their R-hat, ESS
#       and tau, and the acceptance lines
# The references come from the research run, not from this package. --reuse checks an existing --outdir without rerunning.
# (1) involves no floating-point arithmetic and is exact anywhere; (2)-(4) are exact on the reference runtime (Apple M4,
# Julia 1.12.6, the Manifest's versions); elsewhere (2)-(3) are checked to RTOL_PORTABLE and (4) is reported, not failed.
# Mutation test: test/calibration/mutation_test_postprocess.sh.

using CSV, DataFrames, SHA, Printf

const ROOT = dirname(dirname(@__DIR__))
const REF = joinpath(ROOT, "test", "reference")
const SCRIPT = joinpath(ROOT, "scripts", "postprocess_calibration.jl")
const POS = filter(a -> !startswith(a, "--"), ARGS)
length(POS) == 1 || error("usage: gate_postprocess.jl <SLR-RFF-BRICK repo> [--outdir=DIR] [--reuse]")
const RESEARCH = abspath(POS[1])
arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const OUTDIR = abspath(arg("--outdir=", joinpath(ROOT, "outputs", "calibration", "gate_postprocess")))
const REUSE = "--reuse" in ARGS
const SEEDS = (2026, 2027, 2028, 2029)
const N_ITER = 2_000_000
const REFERENCE_RUNTIME = Sys.isapple() && Sys.ARCH === :aarch64 && VERSION == v"1.12.6"
const RTOL_PORTABLE = 1e-10
## Repeated from tools/build_test_reference.py, which extracted the reference lines with them.
const LOG_PATTERNS = [r"^  seed\d+: sd_gis median ", r"^  (seed\d+|POOLED) +-?\d", r"^  SLR@\d+ ", r"<-- check$",
                      r"^  (\d+ params NOT converged\.|all params converged)", r"^\d+ chains × \d+ draws\. Convergence"]

ok = true
fail(msg) = (global ok = false; println("[GATE] FAIL ", msg))

## the inputs are the production chains
sources = Dict(split(l, ",")[2] => split(l, ",")[4] for l in readlines(joinpath(REF, "SOURCES.csv"))[2:end])
for s in SEEDS
    rel = "outputs/mcmc/chain_L27_seed$(s)_n$(N_ITER).csv"
    f = joinpath(RESEARCH, rel)
    isfile(f) || error("missing $f")
    want = sources[rel]
    got = open(io -> bytes2hex(sha256(io)), f)
    got == want || error("$rel: sha256 $got is not the production chain's ($want)")
end
println("[GATE] inputs: the four production chains (sha256 == test/reference/SOURCES.csv)")

## run
mkpath(OUTDIR)
const LOG = joinpath(OUTDIR, "run.log")
if !REUSE
    cmd = `$(Base.julia_cmd()) --project=$ROOT --threads=1 $SCRIPT --chains=$(joinpath(RESEARCH, "outputs", "mcmc")) --outdir=$OUTDIR --accept-slr`
    t = @elapsed run(pipeline(cmd; stdout=LOG, stderr=LOG))
    @printf("[GATE] ran scripts/postprocess_calibration.jl in %.0f s (log %s)\n", t, LOG)
end

## (1) the posterior, byte for byte
post, shipped = joinpath(OUTDIR, "ladrillo_L27_posterior.csv"), joinpath(ROOT, "data", "posterior", "ladrillo_L27_posterior.csv")
if !isfile(post)
    fail("(1) the run wrote no ladrillo_L27_posterior.csv")
elseif read(post) == read(shipped)
    println("[GATE] (1) posterior == data/posterior/ladrillo_L27_posterior.csv, byte for byte ($(filesize(post)) bytes)")
else
    a, b = readlines(post), readlines(shipped)
    i = something(findfirst(k -> a[k] != b[k], 1:min(length(a), length(b))), min(length(a), length(b)) + 1)
    fail("(1) posterior differs from the shipped file: $(length(a)) vs $(length(b)) lines, first difference at line $i")
end

## (2) the certificate
cert, rcert = joinpath(OUTDIR, "slr_convergence_L27.csv"), joinpath(REF, "research_slr_convergence_L27.csv")
if read(cert) == read(rcert)
    println("[GATE] (2) certificate == research slr_convergence_L27.csv, byte for byte")
else
    a, b = CSV.read(cert, DataFrame), CSV.read(rcert, DataFrame)
    num = [c for c in names(b) if eltype(b[!, c]) <: AbstractFloat]
    same_shape = names(a) == names(b) && nrow(a) == nrow(b)
    close_enough = same_shape && all(isapprox(a[!, c], b[!, c]; rtol=RTOL_PORTABLE) for c in num) &&
                   all(a[!, c] == b[!, c] for c in setdiff(names(b), num))
    if !REFERENCE_RUNTIME && close_enough
        println("[GATE] (2) certificate == research to rtol $RTOL_PORTABLE (not the reference runtime)")
    else
        fail("(2) certificate differs from the research file:\n" * read(cert, String) * "research:\n" * read(rcert, String))
    end
end

## (3) the pooled covariance
cv, rcv = joinpath(OUTDIR, "adapted_cov_L27.csv"), joinpath(REF, "research_adapted_cov_L27.csv")
if !isfile(cv)
    fail("(3) the run wrote no adapted_cov_L27.csv")
else
    a, b = CSV.read(cv, DataFrame), CSV.read(rcv, DataFrame)
    names(a) == String.(split(readline(shipped), ",")) || fail("(3) covariance header is not the posterior's header")
    A, B = Matrix{Float64}(a), Matrix{Float64}(b)
    if size(A) != size(B)
        fail("(3) covariance is $(size(A)), research $(size(B))")
    elseif A == B
        println("[GATE] (3) covariance == research adapted_cov_L27.csv at all $(length(A)) entries")
    elseif !REFERENCE_RUNTIME && isapprox(A, B; rtol=RTOL_PORTABLE)
        println("[GATE] (3) covariance == research to rtol $RTOL_PORTABLE (not the reference runtime)")
    else
        fail("(3) covariance differs at $(count(A .!= B)) of $(length(A)) entries; max |rel diff| " *
             "$(maximum(abs.(A .- B) ./ max.(abs.(B), floatmin())))")
    end
end

## (4) the log lines
got = [l for l in readlines(LOG) if any(p -> occursin(p, l), LOG_PATTERNS)]
want = readlines(joinpath(REF, "research_postprocess_L27_log_lines.txt"))
if got == want
    println("[GATE] (4) the run's $(length(got)) diagnostic log lines == the research log's")
else
    msg = "(4) log lines differ ($(length(got)) vs research $(length(want))):\n" *
          join(["  - $w\n  + $g" for (g, w) in zip(got, want) if g != w], "\n") *
          (length(got) == length(want) ? "" : "\n  run:\n" * join(got, "\n"))
    REFERENCE_RUNTIME ? fail(msg) : println("[GATE] (4) REPORTED, not failed (not the reference runtime): ", msg)
end

println(ok ? "[GATE] PASS" : "[GATE] FAIL")
exit(ok ? 0 : 1)
