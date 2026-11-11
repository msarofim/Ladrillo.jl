#!/bin/bash
# Produce every product the identity gates 2-4 compare, with this package's DEFAULTS (v0.2), in three streams.
#   test/identity/run_products.sh            then   test/identity/compare_all.sh <SLR-RFF-BRICK outputs dir>
# Writes outputs/projections/{joint,brick20}_*_<scenario>.csv, fixed_components_ssp126_ssp245_ssp585.csv and
# outputs/hindcast/. Log: outputs/log_run_products.txt.
set -uo pipefail
cd "$(dirname "$0")/../.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
SCEN="ssp126 ssp245 ssp585 vvVL vvL vvLN vvML vvM vvHL vvH"
LOG=outputs/log_run_products.txt; mkdir -p outputs; : > "$LOG"
J="julia --project=."
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
## per-stream log as an argument (macOS /bin/bash is 3.2: no $BASHPID)
step(){ local nm="$1" lg="$2"; shift 2; say "START  $nm"; if "$@" >> "$lg" 2>&1; then say "OK     $nm"; else say "FAILED $nm"; fi; }
s_joint(){ local lg=outputs/log_run_products_joint.txt; : > "$lg"
           for s in $SCEN; do step "joint $s" "$lg" $J scripts/project_joint.jl --scenario=$s; done; }
s_brick(){ local lg=outputs/log_run_products_brick20.txt; : > "$lg"
           for s in $SCEN; do step "brick20 $s" "$lg" $J scripts/brick20_joint.jl --scenario=$s; done; }
s_rest(){  local lg=outputs/log_run_products_rest.txt; : > "$lg"
           step "fixed panel" "$lg" $J scripts/project_fixed.jl
           step "posterior predictive" "$lg" $J scripts/posterior_predictive.jl
           step "BRICK 2.0 posterior predictive" "$lg" $J scripts/posterior_predictive_brick20.jl
           step "IC residuals" "$lg" $J scripts/ic_hindcast_residuals.jl; }
say "products | Ladrillo $(git rev-parse --short HEAD) | load $(uptime | sed 's/.*load averages*: //')"
s_joint & P1=$!; s_brick & P2=$!; s_rest & P3=$!
wait $P1 $P2 $P3
say "DONE: $(grep -c '^.*OK ' "$LOG") ok, $(grep -c 'FAILED' "$LOG") failed"
