#!/bin/bash
## Mutation test of the post-processing gate (test/calibration/gate_postprocess.jl). Each mutant is applied to a scratch
## COPY of the package (the working tree is never edited) and the gate must FAIL. A gate that passes is not a gate that
## works. One mutant per part of the gate, each a slip that leaves the other parts unchanged:
##   P1  (1) subsample takes the LAST draw of each 400-draw block instead of the first (pool[step:step:end])
##   P2  (2) certificate thins from the second draw of each chain's second half (step:step instead of 1:step)
##   P3  (3) covariance jitter 1e-10 -> 1e-9
##   P4  (4) parameter ESS maxlag cap 200,000 -> 100,000
## NO-POWER control, expected to SURVIVE:
##   P5  noise-mode threshold 0.10 -> 0.20 cm: every chain's median sd_gis is ~0.025 cm, so no threshold above that
##       changes anything the gate sees
##   bash test/calibration/mutation_test_postprocess.sh <SLR-RFF-BRICK repo>     (about 25 min per mutant)
## Writes test/calibration/mutation_result_postprocess.txt.
set -uo pipefail
RESEARCH="$(cd "$1" && pwd)"
RESULT="$(cd "$(dirname "$0")" && pwd)/mutation_result_postprocess.txt"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
declare -a NAMES=(P1 P2 P3 P4 P5)
declare -a FILES=(src/postprocess.jl src/postprocess.jl src/postprocess.jl src/postprocess.jl src/postprocess.jl)
declare -a FROM=('draws = pool[1:step:end, :]'
                 'rows = collect(1:step:nrow(ch))'
                 'const COV_JITTER = 1e-10'
                 'const ESS_MAXLAG = 200_000'
                 'const NOISE_SD_GIS_MAX = 0.10 ')
declare -a TO=('draws = pool[step:step:end, :]'
               'rows = collect(step:step:nrow(ch))'
               'const COV_JITTER = 1e-9'
               'const ESS_MAXLAG = 100_000'
               'const NOISE_SD_GIS_MAX = 0.20 ')
echo "post-processing gate mutation test, $(date '+%Y-%m-%d %H:%M'), julia $(julia --version | cut -d' ' -f3)" > "$RESULT"
for i in "${!NAMES[@]}"; do
  rm -rf "$WORK/pkg"; mkdir -p "$WORK/pkg/test"
  cp -R "$ROOT/src" "$ROOT/scripts" "$ROOT/data" "$ROOT/Project.toml" "$ROOT/Manifest.toml" "$WORK/pkg/"
  cp -R "$ROOT/test/calibration" "$ROOT/test/reference" "$WORK/pkg/test/"
  F="${FILES[$i]}"
  python3 - "$WORK/pkg/$F" "${FROM[$i]}" "${TO[$i]}" <<'EOF'
import sys
p, a, b = sys.argv[1:]
s = open(p).read()
n = s.count(a)
if n != 1:
    sys.exit("pattern found %d times" % n)
open(p, "w").write(s.replace(a, b))
EOF
  if cmp -s "$WORK/pkg/$F" "$ROOT/$F"; then echo "${NAMES[$i]}: mutation did not apply -- mutant NOT tested" | tee -a "$RESULT"; continue; fi
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    julia --project="$WORK/pkg" --threads=1 "$WORK/pkg/test/calibration/gate_postprocess.jl" "$RESEARCH" \
    > "$WORK/out_${NAMES[$i]}.txt" 2>&1
  rc=$?
  why="$(grep -E '^\[GATE\] FAIL \(' "$WORK/out_${NAMES[$i]}.txt" | cut -c1-110 | tr '\n' ' ')"
  if [ $rc -eq 0 ]; then echo "${NAMES[$i]}: gate PASSED -- mutant SURVIVED" | tee -a "$RESULT"
  elif grep -q '^\[GATE\] FAIL$' "$WORK/out_${NAMES[$i]}.txt"; then echo "${NAMES[$i]}: gate FAILED -- mutant KILLED  ($why)" | tee -a "$RESULT"
  else echo "${NAMES[$i]}: gate did not run (exit $rc):" | tee -a "$RESULT"; tail -5 "$WORK/out_${NAMES[$i]}.txt" | tee -a "$RESULT"; fi
done
