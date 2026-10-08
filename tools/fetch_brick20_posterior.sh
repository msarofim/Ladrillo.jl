#!/bin/bash
## Fetch the BRICK 2.0 comparison posterior (MimiBRICK v2.0.0, post-PR#93 posterior subsample, 10,000 draws) into
## data/brick20/. It is NOT redistributed in this repository (Marcus 2026-10-07); scripts/brick20_joint.jl needs it.
## The checksum is the file every published BRICK 2.0 number in the v1.0 paper was computed from.
set -euo pipefail
cd "$(dirname "$0")/.."
##
## ⚠ NOT the public MimiBRICK v2.0.0 Zenodo file (checked 2026-10-08). Zenodo 10.5281/zenodo.20592337 (2026-06-08)
## ships a parameters_subsample_brick.csv (7,011,352 B, md5 6487623e..., sha256 fecabef3...) that is a DIFFERENT
## calibration: 0 of 350,000 cells match this file, antarctic_precip0 is stored in log space, and thermal_alpha's median
## is 0.054 against 0.158 here. On this package's FaIR forcing it gives a thermal-expansion hindcast RMSE of 1.97 cm
## (this file: 0.52) and an SSP2-4.5 2100 TE of 6.6 cm (19.4). Earlier versions of that record (17966962, 6626335,
## 6461560) do not match either. This file is Tony Wong's post-PR#93 subsample of 2026-05-22. It has no public URL
## yet.
URL="${BRICK20_POSTERIOR_URL:-}"          # TODO: a public URL for THIS file (sha256 below); none exists yet
SHA256=78247db8babf611151c58149956e3ddd5f6f5f94ce460b0a6bf50c07e9158b35
DEST=data/brick20/parameters_subsample_brick.csv
[[ -n "$URL" ]] || { echo "no public URL yet: set BRICK20_POSTERIOR_URL, or place the file at $DEST"; exit 1; }
mkdir -p data/brick20
curl -fsSL "$URL" -o "$DEST"
echo "$SHA256  $DEST" | shasum -a 256 -c -
