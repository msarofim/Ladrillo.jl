#!/bin/bash
## FACTS on the three SSP controls, TO 2300 (2026-09-11), for the L24 deliverable's FIG 8 (the
## SSP comparison at 2300 had no FACTS bars -- only the van Vuuren set had been extended, see
## run_vv_facts_2300.sh). Same pattern: separate *2300 keys, the 2150 controls untouched so the
## shipped 2100/2150 comparison column stays like-for-like. No emulandice on these keys
## (2100-capped, and per-SSP-trained from a label these keys do not carry).
set -uo pipefail
cd "$(dirname "$0")"
REPO="$(pwd)"
LOG=ssp_facts_2300.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
for e in global.shared.ssp1262300.n200 global.shared.ssp2452300.n200 global.shared.ssp5852300.n200; do
  say "START $e"
  if docker run --rm --volume="$REPO":/opt/facts -w /opt/facts facts \
       /bin/bash -lc "python3 runFACTS.py experiments/$e" >> "$LOG" 2>&1
  then say "OK    $e"; else say "FAILED $e (rc=$?)"; fi
done
say "DONE"
