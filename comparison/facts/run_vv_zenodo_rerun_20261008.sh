#!/bin/bash
## FACTS's seven van Vuuren experiments (+ their *2300 twins, + the vvHL2300 crate0 receipt arm) RE-RUN ON THE
## PUBLISHED EMISSIONS (Marcus 2026-10-08: "go with b"; SLR-RFF-BRICK CHANGELOG 2026-10-08j, handoff 10-08c).
##
## WHY. The vv FaIR cubes were built from the `harmonized` emissions: the IIASA prerelease to 2100 with the published
## extension ADDED at the prerelease 2100 level, which for Medium-to-Low left the CO2 tail 10.5 GtCO2/yr low for 200
## years (0.73 K too cold at 2300). SLR-RFF-BRICK/run_vv_zenodo_rerun_20261008.sh `cubes` rebuilt the cubes on variant
## b (CMIP7 1.6.0 history + the published 2024-2300 scenario); FACTS reads those cubes through
## build_shared_climate_nc.py, so every vv experiment moves.
##
## WHAT: 14 experiments = 7 base keys (to 2150, emulandice on) + 7 *2300 twins, plus
## global.shared.vvHL2300.n200.crate0 (a hand copy of vvHL2300 with GrIS1f crateyear_end: 0; the receipt for the GMD
## draft's FACTS "112 cm", python/diag_facts_gis_extrap_receipt.py). crate0 was NOT in the 10-07 rerun and was still
## on the Smith basis (input dated 10-01); it gets the rebuilt vvHL2300 input/ and is re-run here.
## NOT touched: the SSP controls; every pulse experiment (the pulse arc's downstream waits, Marcus 10-08).
## QUARANTINE: input/ + output/ of all 15 experiments and the shipped facts_components_shared_n200.csv + the receipt
## CSV move to $Q/ before the rebuild.
## ⚠ build_shared_climate_nc.py writes Source="FaIR-2.2.4-calib1.6.0-CMIP7-injected" whatever the emissions variant,
## so the new NetCDFs carry the SAME Source attribute as the old ones: provenance is this script, its log and the
## quarantine README, not the attribute.
## GATES: [VV-ZENODO-BASIS] before anything moves; [CONFIG]; the extractor's [DRIVER] and [EXTEND]; [SSP-ROWS].
## Needs Colima. ~95 min (10-07: 93 min).
set -uo pipefail
cd "$(dirname "$0")"
REPO="$(pwd)"
SLR=~/Documents/2026/CodeProjects/SLR-RFF-BRICK
SHIPPED=$SLR/outputs/facts_components_shared_n200.csv
RECEIPT=$SLR/outputs/diag_facts_gis_extrap_receipt_vvHL2300.csv
QC=$SLR/outputs/quarantine/20261008_vv_harmonized_tail_cubes/data_observations
ARMB=$SLR/outputs/vvz_regression_20261008/arm_b_cubes
Q=quarantine/20261008_vv_harmonized_tail
LOG=vv_zenodo_rerun_20261008.log; : > "$LOG"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }
BASE="vvVL vvLN vvL vvML vvM vvHL vvH"
KEYS="$BASE $(for k in $BASE; do printf '%s2300 ' $k; done)"
CRATE0=global.shared.vvHL2300.n200.crate0
source ~/climate-env/bin/activate
say "FACTS vv zenodo rerun | facts $(git rev-parse --short HEAD) | SLR $(git -C $SLR rev-parse --short HEAD) | load $(uptime | sed 's/.*load averages*: //')"

python $SLR/python/gate_vv_zenodo_basis.py "$QC" "$ARMB" >> "$LOG" 2>&1 || { say "[VV-ZENODO-BASIS] FAILED -- STOP"; exit 1; }
say "[VV-ZENODO-BASIS] PASS"
docker info > /dev/null 2>&1 || { say "docker is not up (colima start) -- STOP"; exit 1; }

[[ -e "$Q" ]] && { say "$Q exists -- refusing to overwrite a quarantine. STOP"; exit 1; }
mkdir -p "$Q"
cp -p "$SHIPPED" "$Q/"; [[ -f "$RECEIPT" ]] && cp -p "$RECEIPT" "$Q/"
for e in $(for k in $KEYS; do echo global.shared.$k.n200; done) $CRATE0; do
  mkdir -p "$Q/experiments/$e"; mv "experiments/$e/input" "experiments/$e/output" "$Q/experiments/$e/"
done
say "quarantined input/ + output/ of $(( $(echo $KEYS | wc -w) + 1 )) experiments -> $Q"

for k in $KEYS; do
  python build_shared_climate_nc.py --scen=$k >> "$LOG" 2>&1 || { say "FAILED climate build $k -- STOP"; exit 1; }
done
python build_shared_configs.py >> "$LOG" 2>&1 || { say "FAILED build_shared_configs -- STOP"; exit 1; }
git diff --quiet -- experiments/ || { say "[CONFIG] a tracked config changed -- STOP (git diff experiments/)"; exit 1; }
say "[CONFIG] all tracked configs byte-identical"
cp -Rp experiments/global.shared.vvHL2300.n200/input experiments/$CRATE0/input
say "crate0: input/ copied from the rebuilt global.shared.vvHL2300.n200"

FAILS=0
for e in $(for k in $KEYS; do echo global.shared.$k.n200; done) $CRATE0; do
  say "START $e"
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
python $SLR/python/diag_facts_gis_extrap_receipt.py >> "$LOG" 2>&1 && say "receipt (112 cm) re-run: OK" || { say "FAILED receipt"; exit 1; }
say "DONE. Next: SLR-RFF-BRICK/run_vv_zenodo_rerun_20261008.sh arms (after MAGICC-SLR vv)"
