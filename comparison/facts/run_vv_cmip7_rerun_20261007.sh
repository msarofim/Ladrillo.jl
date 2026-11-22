#!/bin/bash
## FACTS's seven van Vuuren experiments (+ their *2300 twins) RE-RUN ON THE CMIP7 FORCING BASIS (Marcus 2026-10-07,
## SLR-RFF-BRICK notes/handoff_2026-10-07_clean_l27_repo.md §2.1: "Rerun all on CMIP7", FACTS included).
##
## WHY. The inputs were built 09-03 from the Smith-history vv cubes: rebuilding
## experiments/global.shared.vvM.n200/input/shared_vvM_gsat.nc with build_shared_climate_nc.py's own splice and
## subsample reproduces it from SLR-RFF-BRICK/outputs/quarantine/20261002_vv_smith_history_basis/ at max|d| = 0.0 K
## and misses the live CMIP7 cubes by up to 0.12 K. The Ladrillo and BRICK 2.0 vv arms are re-run on CMIP7 by
## SLR-RFF-BRICK/run_vv_cmip7_rerun_20261007.sh; without this FACTS column the comparison would be mixed-basis.
##
## WHAT: 14 experiments = 7 base keys (to 2150, emulandice on) + 7 *2300 twins. The SSP controls and every pulse
## experiment are NOT touched (SSP cubes were already CMIP7; the pulse arms belong to the pulse comparison).
## QUARANTINE: each experiment's input/ and output/ move to $Q/<exp>/ before the rebuild (configs stay, and must
## come back byte-identical from build_shared_configs.py since no NetCDF attribute changes); the shipped
## SLR-RFF-BRICK/outputs/facts_components_shared_n200.csv is copied there too.
## GATES after extraction: [SSP-ROWS] every non-vv row of the shipped CSV must be byte-identical to before;
## the extractor's own [DRIVER] and [EXTEND] gates still apply.
## Needs Colima (profile `default`: 8 CPU, 8 GiB, vz, virtiofs). ~45 min.
set -uo pipefail
cd "$(dirname "$0")"
REPO="$(pwd)"
SLR=~/Documents/2026/CodeProjects/SLR-RFF-BRICK
SHIPPED=$SLR/outputs/facts_components_shared_n200.csv
Q=quarantine/20261007_vv_smith_history
LOG=vv_cmip7_rerun_20261007.log; : > "$LOG"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }
BASE="vvVL vvLN vvL vvML vvM vvHL vvH"
KEYS="$BASE $(for k in $BASE; do printf '%s2300 ' $k; done)"
source ~/climate-env/bin/activate
say "FACTS vv CMIP7 rerun | facts $(git rev-parse --short HEAD) | SLR $(git -C $SLR rev-parse --short HEAD)"

python $SLR/python/gate_vv_cmip7_basis.py >> "$LOG" 2>&1 || { say "[VV-BASIS] FAILED -- STOP"; exit 1; }
say "[VV-BASIS] PASS"
docker info > /dev/null 2>&1 || { say "docker is not up (colima start) -- STOP"; exit 1; }

[[ -e "$Q" ]] && { say "$Q exists -- refusing to overwrite a quarantine. STOP"; exit 1; }
mkdir -p "$Q"
cp -p "$SHIPPED" "$Q/"
for k in $KEYS; do
  e=experiments/global.shared.$k.n200
  mkdir -p "$Q/$e"; mv "$e/input" "$e/output" "$Q/$e/"
done
say "quarantined input/ + output/ of $(echo $KEYS | wc -w | tr -d ' ') experiments -> $Q"

for k in $KEYS; do
  python build_shared_climate_nc.py --scen=$k >> "$LOG" 2>&1 || { say "FAILED climate build $k -- STOP"; exit 1; }
done
python build_shared_configs.py >> "$LOG" 2>&1 || { say "FAILED build_shared_configs -- STOP"; exit 1; }
git diff --quiet -- experiments/ || { say "[CONFIG] a tracked config changed -- STOP (git diff experiments/)"; exit 1; }
say "[CONFIG] all tracked configs byte-identical"

FAILS=0
for k in $KEYS; do
  e=global.shared.$k.n200; say "START $e"
  if docker run --rm --volume="$REPO":/opt/facts -w /opt/facts facts \
       /bin/bash -lc "python3 runFACTS.py experiments/$e" >> "$LOG" 2>&1
  then say "OK    $e"; else say "FAILED $e (rc=$?)"; FAILS=$((FAILS+1)); fi
done
[[ $FAILS -eq 0 ]] || { say "$FAILS experiment(s) failed -- not extracting. STOP"; exit 1; }

python extract_facts_shared_components.py >> "$LOG" 2>&1 || { say "FAILED extraction -- STOP"; exit 1; }
python - "$Q/$(basename $SHIPPED)" "$SHIPPED" >> "$LOG" 2>&1 <<'EOF' || { say "[SSP-ROWS] FAILED"; exit 1; }
import sys, pandas as pd
old, new = (pd.read_csv(f, dtype=str) for f in sys.argv[1:3])
col = next(c for c in ("scenario", "key", "ssp") if c in old.columns)
o, n = (d[~d[col].str.startswith("vv")].reset_index(drop=True) for d in (old, new))
ok = o.equals(n)
print("[SSP-ROWS] %d non-vv rows old vs new: %s" % (len(o), "IDENTICAL" if ok else "DIFFER"))
sys.exit(0 if ok else 1)
EOF
say "[SSP-ROWS] PASS"
say "DONE. Next: SLR-RFF-BRICK/run_vv_cmip7_rerun_20261007.sh downstream"
