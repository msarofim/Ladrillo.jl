# Calibration gate: does src/calibration.jl define the production L27 objective, and does scripts/calibrate.jl's sampler
# reproduce the production chain?
#
#   julia --project=. --threads=1 test/calibration/gate_chain_prefix.jl [n_rows]      (default 20000, about 2.5 min)
#
# Exit 0 = PASS. The references are the research production chains (SLR-RFF-BRICK outputs/mcmc/chain_L27_seed*_n2000000.csv,
# the run behind the shipped posterior), copied as text by tools/build_test_reference.py; their sha256 are in
# test/reference/SOURCES.csv. They come from the research calibrator, not from this package, so a pass certifies the
# objective rather than re-baselining it.
#   (1) logposterior at each of the four start rows == the first row of that seed's production chain (log_post field)
#   (2) seed 2026 for n_rows iterations: the log_post of EVERY row == production, and the full row (all 50 parameters)
#       == production at the reference rows (1, 459-461, every 500th)
# WHY 20,000 ROWS. The production chain's first accepted move is row 460: a shorter window (the research gate used 300)
# compares nothing but the start and a run of rejections. The first 20,000 rows hold 4,603 accepted moves.
# EXACT ONLY ON THE REFERENCE RUNTIME (Apple M4, Julia 1.12.6, the Manifest's package versions, one BLAS thread). Elsewhere
# a last-bit difference can flip one accept decision and the chains part; (1) is then checked to RTOL_PORTABLE and (2) is
# reported, not failed.
# Mutation test: test/calibration/mutation_test.sh.

using Ladrillo, Printf
import LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

const C = Ladrillo.Calibration
const REF = joinpath(dirname(@__DIR__), "reference")
const N_ROWS = isempty(ARGS) ? 20000 : parse(Int, ARGS[1])
const SEED = 2026
const REFERENCE_RUNTIME = Sys.isapple() && Sys.ARCH === :aarch64 && VERSION == v"1.12.6" &&
                          Base.JLOptions().check_bounds == 0
const RTOL_PORTABLE = 1e-12
"""The research chains name the glacier reservoirs SLOWP and FAST; the package names them SLOWG and FASTG."""
pkgname(s) = replace(replace(s, r"_SLOWP$" => "_SLOWG"), r"_FAST$" => "_FASTG")

function read_ref(name)
    lines = readlines(joinpath(REF, name))
    return String.(split(lines[1], ",")), [String.(split(l, ",")) for l in lines[2:end]]
end

ok = true
p = C.problem()
## (1) the start rows of the four production chains
for seed in C.SEEDS_CANON
    hdr, rows = read_ref("research_chain_L27_seed$(seed)_row1.csv")
    pkgname.(hdr[1:end-2]) == C.param_names(p) || error("seed $seed: reference header does not match the parameter names")
    θref = parse.(Float64, rows[1][1:end-2])
    θ0 = C.start_point(p, seed)
    θ0 == θref || (global ok = false; println("[GATE] FAIL seed $seed: start row differs from the production chain's row 1"))
    lp, lpref = C.logposterior(p, θ0), parse(Float64, rows[1][end-1])
    pass = REFERENCE_RUNTIME ? lp == lpref : isapprox(lp, lpref; rtol=RTOL_PORTABLE)
    pass || (global ok = false)
    @printf("[GATE] (1) seed %d  logposterior(start) %s  production %s  %s\n", seed, repr(lp), repr(lpref),
            pass ? "PASS" : "FAIL")
end

## (2) the first N_ROWS of the seed-2026 chain
_, lpr = read_ref("research_chain_L27_seed2026_prefix_logpost.csv")
N_ROWS <= length(lpr) || error("the reference holds $(length(lpr)) rows; asked for $N_ROWS")
t = @elapsed r = C.run_chain(p, N_ROWS, SEED; show_progress=false)
lpref = [parse(Float64, x[2]) for x in lpr[1:N_ROWS]]
nmoves = count(i -> r.logpost[i] != r.logpost[i-1], 2:N_ROWS)
bad = findfirst(i -> r.logpost[i] != lpref[i], 1:N_ROWS)
@printf("[GATE] (2) seed %d, %d iterations in %.0f s, acceptance %.4f, %d moves (production prefix: %d)\n", SEED, N_ROWS, t,
        r.accept, nmoves, count(i -> lpref[i] != lpref[i-1], 2:N_ROWS))
if bad === nothing
    println("[GATE] (2) log_post == production at all $N_ROWS rows")
else
    println("[GATE] (2) log_post first differs at row $bad: $(repr(r.logpost[bad])) vs production $(repr(lpref[bad]))")
end
hdr, rows = read_ref("research_chain_L27_seed2026_prefix_rows.csv")
pkgname.(hdr[1:end-2]) == C.param_names(p) || error("prefix rows: reference header does not match the parameter names")
nrow_checked, nrow_bad = 0, 0
## which production rows the reference holds: recorded by tools/build_test_reference.py in SOURCES.csv
const CHAIN_ROWS = let l = only(filter(l -> startswith(l, "research_chain_L27_seed2026_prefix_rows.csv,"),
                                       readlines(joinpath(REF, "SOURCES.csv"))))
    parse.(Int, split(split(l, ",")[end]))
end
length(CHAIN_ROWS) == length(rows) || error("prefix rows: SOURCES.csv lists $(length(CHAIN_ROWS)) rows, the file has $(length(rows))")
for (i, row) in zip(CHAIN_ROWS, rows)
    parse(Float64, row[end-1]) == lpref[i] || error("prefix rows: row $i's log_post disagrees with the log_post reference")
    i <= N_ROWS || continue
    global nrow_checked += 1
    r.chain[i, :] == parse.(Float64, row[1:end-2]) || (global nrow_bad += 1)
end
pass2 = bad === nothing && nrow_bad == 0
println("[GATE] (2) full rows == production at $(nrow_checked - nrow_bad) of $nrow_checked reference rows")
if REFERENCE_RUNTIME
    ok &= pass2
    println(ok ? "[GATE] PASS: the package calibrator reproduces the production L27 objective and the first $N_ROWS rows of its chain" :
                 "[GATE] FAIL")
else
    println("[GATE] not the reference runtime: (2) is reported, not gated" * (pass2 ? " (it is identical here too)" : ""))
    println(ok ? "[GATE] PASS (portable): logposterior(start) agrees to $RTOL_PORTABLE for all four seeds" : "[GATE] FAIL")
end
exit(ok ? 0 : 1)
