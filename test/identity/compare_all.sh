#!/bin/bash
# Identity gates 2-4 against the research outputs, after test/identity/run_products.sh.
#   test/identity/compare_all.sh <SLR-RFF-BRICK outputs dir>
# Gate 1 (kernel) is dump_old.jl + dump_new.jl + compare.jl; gate 5 (Table 5) is compare_ic_table.py, needed only when
# the residuals move (gate 4 says whether they did).
set -uo pipefail
RES="${1:?usage: $0 <SLR-RFF-BRICK outputs dir>}"
cd "$(dirname "$0")/../.."
SCEN="ssp126 ssp245 ssp585 vvVL vvL vvLN vvML vvM vvHL vvH"
J="julia --project=."
fails=0
run(){ echo "== $*"; "$@" || fails=$((fails+1)); }
run $J test/identity/compare_products.jl "$RES" ladrillo $SCEN
run $J test/identity/compare_products.jl "$RES" brick20 $SCEN
run $J test/identity/compare_panel.jl "$RES"
run $J test/identity/compare_postpred.jl "$RES"
run $J test/identity/compare_hindcast.jl "$RES"
echo "[GATES 2-4] $([[ $fails == 0 ]] && echo PASS || echo "FAIL ($fails gate(s))")"
exit $fails
