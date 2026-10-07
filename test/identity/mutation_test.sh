#!/bin/bash
## Mutation test of the v1.0 identity gate: each mutant must make test/identity/compare.jl FAIL on a case whose
## reference (old_<case>.jls) already exists in $LADRILLO_IDENTITY_DIR. A gate that passes is not a gate that works.
##   test/identity/mutation_test.sh <case>      (K3_ssp245harm_central is the cheap one; K4 exercises the paleo block)
## Mutants (each applied to a scratch COPY of the package; the working tree is never edited):
##   M1  active-basin share computed as 4.62/7.35 -- the same number on paper, possibly not in the last bit
##   M2  the Greenland headroom summed in a different order -- float associativity
##   M3  paleo rows attached in one block instead of the case's blocks (only bites on K4)
##   M4  the glacier floor 1e-12 -> 1e-11 (a no-power mutant: expected to SURVIVE if the floor never binds)
set -uo pipefail
CASE="${1:?case}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
: "${LADRILLO_IDENTITY_DIR:?set LADRILLO_IDENTITY_DIR to where old_<case>.jls lives}"
[[ -f "$LADRILLO_IDENTITY_DIR/old_$CASE.jls" ]] || { echo "no reference old_$CASE.jls"; exit 1; }
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
declare -a NAMES=(M1 M2 M3 M4)
declare -a FILES=(src/components/greenland.jl src/components/greenland.jl src/posterior.jl src/components/glaciers.jl)
declare -a FROM=('GIS_SECTOR_VOL_M.south \/ sum(GIS_SECTOR_VOL_M) + GIS_SECTOR_VOL_M.mid \/ sum(GIS_SECTOR_VOL_M)'
                 'head = max(p.gis_v0 - (fa + fh + sa + sh), 0.0)'
                 'attach_paleo!(df; block = paleo_block === nothing ? nrow(df) : paleo_block)'
                 'frac_left = max(1 - S \/ a, 1e-12)')
declare -a TO=('4.62 \/ 7.35'
               'head = max(p.gis_v0 - (fa + sa + fh + sh), 0.0)'
               'attach_paleo!(df)'
               'frac_left = max(1 - S \/ a, 1e-11)')
for i in "${!NAMES[@]}"; do
  rm -rf "$WORK/pkg"; cp -R "$ROOT" "$WORK/pkg"; rm -rf "$WORK/pkg/.git"
  f="$WORK/pkg/${FILES[$i]}"
  sed -i '' "s/${FROM[$i]}/${TO[$i]}/" "$f"
  if cmp -s "$f" "$ROOT/${FILES[$i]}"; then echo "${NAMES[$i]}: sed did not apply -- mutant NOT tested"; continue; fi
  OUT="$WORK/out_${NAMES[$i]}"; mkdir -p "$OUT"; cp "$LADRILLO_IDENTITY_DIR/old_$CASE.jls" "$OUT/"
  LADRILLO_IDENTITY_DIR="$OUT" julia --project="$WORK/pkg" "$WORK/pkg/test/identity/dump_new.jl" "$CASE" > /dev/null 2>&1 \
    || { echo "${NAMES[$i]}: dump failed"; continue; }
  if LADRILLO_IDENTITY_DIR="$OUT" julia --project="$WORK/pkg" "$WORK/pkg/test/identity/compare.jl" "$CASE" > "$OUT/cmp.txt" 2>&1
  then echo "${NAMES[$i]}: gate PASSED -- mutant SURVIVED"
  else echo "${NAMES[$i]}: gate FAILED -- mutant KILLED  ($(grep -c FAIL "$OUT/cmp.txt") component rows failed)"; fi
done
