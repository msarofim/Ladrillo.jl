#!/bin/bash
## Mutation test of the calibration gate (test/calibration/gate_chain_prefix.jl). Each mutant is applied to a scratch
## COPY of the package (src/, scripts/, test/calibration/, test/reference/, data/, the project files; the working tree is
## never edited) and the gate must FAIL. A gate that passes is not a gate that works.
##   bash test/calibration/mutation_test.sh [n_rows]      (default 20000; about 2 min per mutant)
## Mutants that leave logposterior at all four start rows UNCHANGED, so only part (2), the chain, can kill them (the
## research repository's 300-row gate, whose window ends before the first accepted move at row 460, misses all three):
##   M1  sampler target acceptance 0.234 -> 0.235
##   M2  AR(1) rho bound 0.99 -> 0.98 (the four start rows' largest rho is 0.976; the seed-2026 chain reaches 0.98 by
##       row 1,000). A cap of this kind went undetected by the research gate in September.
##   M4  the active-basin share targets written as 0.799 / 0.816, equal on paper to 0.592 + 0.207 and 0.554 + 0.262: a
##       last-bit change to a likelihood constant
## Mutant that moves logposterior(start), killed by part (1):
##   M3  GlaMBIE partition sd 0.05 -> 0.0501
## NO-POWER control, expected to SURVIVE:
##   M5  the geometry block's diagonal proposal scale 0.02 -> 0.03: every geometry entry is then overwritten by the
##       L26 covariance seed, so the constant does not reach the sampler
## Writes test/calibration/mutation_result.txt.
set -uo pipefail
NROWS="${1:-20000}"
RESULT="$(cd "$(dirname "$0")" && pwd)/mutation_result.txt"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
F=src/calibration.jl
declare -a NAMES=(M1 M2 M3 M4 M5)
declare -a FROM=('const OPT_ALPHA = 0.234'
                 'const RHO_MAX = 0.99 '
                 'const GLAMBIE_SHARE_SD = 0.05'
                 'const GISB_SHARE_ACTIVE = Tuple(sh.south + sh.mid for sh in GISB_SHARE3)'
                 'const PROP_SCALE, GEO_PROP_SCALE = 0.1, 0.02')
declare -a TO=('const OPT_ALPHA = 0.235'
               'const RHO_MAX = 0.98 '
               'const GLAMBIE_SHARE_SD = 0.0501'
               'const GISB_SHARE_ACTIVE = (0.799, 0.816)'
               'const PROP_SCALE, GEO_PROP_SCALE = 0.1, 0.03')
echo "calibration gate mutation test, $NROWS rows, $(date '+%Y-%m-%d %H:%M'), julia $(julia --version | cut -d' ' -f3)" > "$RESULT"
for i in "${!NAMES[@]}"; do
  rm -rf "$WORK/pkg"; mkdir -p "$WORK/pkg/test"
  cp -R "$ROOT/src" "$ROOT/scripts" "$ROOT/data" "$ROOT/Project.toml" "$ROOT/Manifest.toml" "$WORK/pkg/"
  cp -R "$ROOT/test/calibration" "$ROOT/test/reference" "$WORK/pkg/test/"
  python3 - "$WORK/pkg/$F" "${FROM[$i]}" "${TO[$i]}" <<'EOF'
import sys
p, a, b = sys.argv[1:]
s = open(p).read()
n = s.count(a)
if n != 1:
    sys.exit("pattern found %d times" % n)
open(p, "w").write(s.replace(a, b))
EOF
  if cmp -s "$WORK/pkg/$F" "$ROOT/$F"; then echo "${NAMES[$i]}: mutation did not apply -- mutant NOT tested"; continue; fi
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    julia --project="$WORK/pkg" --threads=1 "$WORK/pkg/test/calibration/gate_chain_prefix.jl" "$NROWS" > "$WORK/out_${NAMES[$i]}.txt" 2>&1
  rc=$?
  why="$(grep -E 'FAIL|first differs' "$WORK/out_${NAMES[$i]}.txt" | head -2 | tr '\n' ' ')"
  if [ $rc -eq 0 ]; then echo "${NAMES[$i]}: gate PASSED -- mutant SURVIVED" | tee -a "$RESULT"
  elif grep -q '^\[GATE\]' "$WORK/out_${NAMES[$i]}.txt"; then echo "${NAMES[$i]}: gate FAILED -- mutant KILLED  ($why)" | tee -a "$RESULT"
  else echo "${NAMES[$i]}: gate did not run (exit $rc):" | tee -a "$RESULT"; tail -5 "$WORK/out_${NAMES[$i]}.txt"; fi
done
