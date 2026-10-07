# Changelog

## 2026-10-07 — package built from the research repository; certified bit-identical

**What this is.** A clean, L27-only `Ladrillo.jl` (Marcus 2026-10-07), built from `SLR-RFF-BRICK` (`ladrillo-dev`,
commit 0c4de0a). The glacier blocks are renamed SLOWP→SLOWG and FAST→FASTG; the rename is label-only
(`data/MANIFEST.csv`, `tools/import_l27_data.py`).

**Structural simplifications**, each covered by the identity gate below:
- Greenland is a dedicated two-basin component. It replaces the three-basin component run with a zero middle share.
- The shipped tap is now structural: a two-stage cascade with the whole-sheet clamp.
- The glacier regrowth divisor is dropped (it was 1.0).
- Every non-L27 variant path is removed.

**Identity gate 1 (kernel)** compares the package with the research kernel on every component × year × draw, using
`==`. All four cases pass: ssp585 on all 10,000 draws with the tap on; ssp245 with no tap; the hindcast under :central
land water; and the vvH joint-climate arm.

**Identity gate 2 (products)** compares the package's joint-climate drivers with the research drivers' outputs. On
vvVL, all 72,000 Ladrillo values and all 36,000 BRICK 2.0 values are identical. The research driver read the raw
2 GB chains; the package reads the 10k subsample. That proves the decision-2 switch moved nothing:
`subsample[1:5:end]` is bit-identical to the joint arm's chain draws.

**Mutation test** (case K4). Killed: a 1-ulp change to the active-basin share, a re-associated glacier sum, and a
paleo assignment in one block. Survived, as expected: two no-power controls.
- **Tried and replaced:** the first reassociation mutant reordered the Greenland headroom sum. It survived because the
  whole-sheet clamp never binds, so it is kept as a no-power control.

**The MimiBRICK build scenario was checked:** building with ssp585 vs ssp245 changes 0 of 200 draws × 6 components ×
451 years. So every product is built on ssp245, which also lets the van Vuuren scenarios run.

**Decisions (Marcus 2026-10-07):**
- **Paleo attachment.** Keep v1.0's rules: the joint arm assigns per chain (`paleo_block = 500`); the SSP panel
  assigns once over its 2,000 rows. Move to a single assignment at v1.1.
  - Measured on ssp585: medians move ≤ 1.3 bootstrap se.
  - The 2100 p95 moves −3.6 cm (total) and −5.0 cm (AIS), about 2.5 se.
  - Tool: `scripts/project_joint.jl --paleo=single`.
- **BRICK 2.0 posterior.** Not redistributed: `tools/fetch_brick20_posterior.sh` fetches it and checks its sha256,
  but the URL is still to come. The file was purged from this repository's local history before any remote existed.
- **The van Vuuren forcing is the CMIP7 basis** (the research repo's 2026-10-02 regeneration). The research repo's vv
  arms and FACTS are being re-run on it.
