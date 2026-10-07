#!/bin/bash
## Fetch the BRICK 2.0 comparison posterior (MimiBRICK v2.0.0, post-PR#93 posterior subsample, 10,000 draws) into
## data/brick20/. It is NOT redistributed in this repository (Marcus 2026-10-07); scripts/brick20_joint.jl needs it.
## The checksum is the file every published BRICK 2.0 number in the v1.0 paper was computed from.
set -euo pipefail
cd "$(dirname "$0")/.."
URL="${BRICK20_POSTERIOR_URL:-}"          # TODO(v1.0 release): the public URL, from Marcus / Tony Wong
SHA256=78247db8babf611151c58149956e3ddd5f6f5f94ce460b0a6bf50c07e9158b35
DEST=data/brick20/parameters_subsample_brick.csv
[[ -n "$URL" ]] || { echo "no public URL yet: set BRICK20_POSTERIOR_URL, or place the file at $DEST"; exit 1; }
mkdir -p data/brick20
curl -fsSL "$URL" -o "$DEST"
echo "$SHA256  $DEST" | shasum -a 256 -c -
