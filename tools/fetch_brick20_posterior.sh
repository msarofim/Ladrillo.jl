#!/bin/bash
## Fetch the BRICK 2.0 comparison posterior into data/brick20/ and check its sha256. It is NOT redistributed in this
## repository (Marcus 2026-10-07); scripts/brick20_joint.jl needs it.
##
## The file is the joint SNEASY-BRICK posterior subsample published with MimiBRICK v2.0.0: Zenodo
## 10.5281/zenodo.20592337 (2026-06-08, CC-BY-4.0), parameters_subsample_sneasybrick.csv, 10,077,052 B, md5
## 2ac199193ae2139f22d253ba557878e7. The arm reads its 35 BRICK columns (Marcus 2026-10-15).
##
## ⚠ The same record also ships parameters_subsample_brick.csv (sha256 fecabef3...). That is the STANDALONE BRICK
## calibration, driven by one SNEASY trajectory, and it is NOT this arm's file: on this package's FaIR forcing its
## thermal expansion is about 3x too small. Through v1.2.0 the arm used Tony Wong's May-22 subsample (sha256
## 78247db8...), the BRICK columns of an earlier joint posterior; the Zenodo joint file gives an SSP2-4.5 total about
## 4 % lower at 2100 and 2300 (SLR-RFF-BRICK CHANGELOG 2026-10-15c).
set -euo pipefail
cd "$(dirname "$0")/.."
URL="${BRICK20_POSTERIOR_URL:-https://zenodo.org/records/20592337/files/parameters_subsample_sneasybrick.csv?download=1}"
SHA256=d74577f603602b222c829a152cb4dbcb60753af34ce4d1fba4cd0dfea52c203c
DEST=data/brick20/parameters_subsample_sneasybrick.csv
mkdir -p data/brick20
curl -fsSL "$URL" -o "$DEST"
echo "$SHA256  $DEST" | shasum -a 256 -c -
