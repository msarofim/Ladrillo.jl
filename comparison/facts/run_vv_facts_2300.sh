#!/bin/bash
## FACTS on the seven van Vuuren markers, TO 2300 (2026-09-03), for the L24 deliverable's Fig 4.
## Same pattern as run_matched_pair_facts_2300.sh -- separate *2300 keys, existing 2150 arms
## untouched so the shipped comparison column stays like-for-like.
set -uo pipefail
cd "$(dirname "$0")"
REPO="$(pwd)"
LOG=vv_facts_2300.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
for e in global.shared.vvVL2300.n200 global.shared.vvLN2300.n200 global.shared.vvL2300.n200 \
         global.shared.vvML2300.n200 global.shared.vvM2300.n200 global.shared.vvHL2300.n200 \
         global.shared.vvH2300.n200; do
  say "START $e"
  if docker run --rm --volume="$REPO":/opt/facts -w /opt/facts facts \
       /bin/bash -lc "python3 runFACTS.py experiments/$e" >> "$LOG" 2>&1
  then say "OK    $e"; else say "FAILED $e (rc=$?)"; fi
done
say "DONE"
