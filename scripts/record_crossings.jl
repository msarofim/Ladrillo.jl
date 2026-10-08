# Which projection draws fire Antarctic fast dynamics within the observed record, per scenario and arm.
#
#   julia --project=. scripts/record_crossings.jl [--outdir=outputs/conditioning]
#
# For the 2,000 projection draws (posterior rows 1:5:10000, v1.1 paleo single assignment) and every shipped scenario,
# the first year fast dynamics is on (src/record_crossing.jl; the model's own trigger, checked against it draw for draw
# by test/conditioning/check_record_crossing.jl) in three climates:
#   onset_record  the calibration's hindcast climate (RECORD_FORCING): what posterior(; drop_record_crossings=true) tests
#   onset_fixed   the scenario's FaIR ensemble mean (the fixed arm)
#   onset_joint   the draw's own FaIR config, spliced to the mean at 2014, paired exactly as scripts/project_joint.jl
#                 pairs it (the joint arm, the paper's reported band)
# 0 = does not fire by END_YEAR. `config` lets a reader check the pairing against a joint-arm draws file.
#
# Writes <outdir>/record_crossings_projection_draws.csv (+ provenance column). Nothing here changes a projection.

using Ladrillo, CSV, DataFrames, Printf, Random

arg(flag, default) = (i = findfirst(a -> startswith(a, flag), ARGS); i === nothing ? default : ARGS[i][length(flag)+1:end])
const OUTDIR = arg("--outdir=", joinpath(dirname(@__DIR__), "outputs", "conditioning"))
const ROWS = 1:5:10000
const PAIR_SEED = 2026                       # scripts/project_joint.jl
const END_YEAR = record_end() + 1            # one year past the record, so a draw first firing just after it shows
const NEVER = 0
const PKG_COMMIT = try readchomp(`git -C $(dirname(@__DIR__)) rev-parse --short HEAD`) catch; "no git" end
const PROVENANCE = "Ladrillo $(pkgversion(Ladrillo)) @ $PKG_COMMIT (posterior L27) | scripts/record_crossings.jl | draws rows $(ROWS) " *
    "(paleo single assignment over 10000, seed $(Ladrillo.PALEO_SEED)) | pairing seed $PAIR_SEED | " *
    "FaIR 2.2.4 (calib 1.6.0) CMIP7-basis forcing, spliced at $(Ladrillo.SPLICE_YEAR) | record climate " *
    "$(Ladrillo.RECORD_FORCING) mean | onset = first year fast dynamics is on, $(Ladrillo.RECORD_START)-$END_YEAR " *
    "($NEVER = not by $END_YEAR) | Antarctic record ends $(record_end()) | julia $(VERSION)"

mkpath(OUTDIR)
post = posterior(rows = ROWS)
const N = nrow(post)
yrs = collect(Ladrillo.RECORD_START:END_YEAR)
enc(y) = y == typemax(Int) ? NEVER : y
onset_record = enc.(record_crossings(post; last_year = END_YEAR))
out = DataFrame[]
for scen in Ladrillo.SCENARIOS
    mg, _ = mean_forcing(scen; years = 1850:2300)
    cy, cfg, G, _ = forcing_cube(scen)
    cy == collect(1850:2300) || error("$scen cube years are not 1850-2300")
    iend = findfirst(==(END_YEAR), cy)
    assign = let rng = MersenneTwister(PAIR_SEED), a = Int[]
        while length(a) < N; append!(a, randperm(rng, length(cfg))); end
        a[1:N]
    end
    gfix = Float64.(mg)[1:iend]
    of = [enc(fastdyn_onset_year(post[k, :], gfix, yrs)) for k in 1:N]
    oj = [enc(fastdyn_onset_year(post[k, :], splice(G[:, assign[k]], Float64.(mg), cy)[1:iend], yrs)) for k in 1:N]
    push!(out, DataFrame(scenario = scen, draw = 1:N, file_row = collect(ROWS), config = cfg[assign],
                         onset_record = onset_record, onset_fixed = of, onset_joint = oj))
    @printf("%-7s fire by %d: record climate %d, fixed %d, joint %d (by %d: %d / %d / %d)\n", scen, record_end(),
            count(y -> 0 < y <= record_end(), onset_record), count(y -> 0 < y <= record_end(), of),
            count(y -> 0 < y <= record_end(), oj), END_YEAR, count(>(0), onset_record), count(>(0), of), count(>(0), oj))
end
df = vcat(out...)
df.provenance .= PROVENANCE
f = joinpath(OUTDIR, "record_crossings_projection_draws.csv")
CSV.write(f, df)
println("wrote ", f)
