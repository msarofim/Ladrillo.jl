#!/usr/bin/env python3
"""
ladrillo_table_a2.py -- Table A2 for the GMD paper: the Antarctic-block combinations the observations identify,
from the principal components written by diag_ais_block_pca.py (run that first). A principal component counts as
"identified" when its posterior variance is below IDENT_MAX of its prior variance; the number of
prior-propagated components (ratio above PROP_MIN) is stated in the caption.

  python python/tables/ladrillo_table_a2.py [--outdir=outputs/tables] [--pca=<outdir>/diag_ais_block_pca.csv]

PORT of SLR-RFF-BRICK python/ladrillo_table_a2.py --tag=L27 --source=subsample. The table code is copied unchanged;
only the I/O differs. Writes <outdir>/ladrillo_table_a2.md (research: outputs/ladrillo_table_a2_L27_sub10k.md).
python/tables/check_tables_identity.py gates it against the research file.
"""
import argparse
import os

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = "python/tables/ladrillo_table_a2.py"
TAG = "L27"                                       # the posterior's calibration
PCA_NAME = "diag_ais_block_pca.csv"               # diag_ais_block_pca.py's OUT_NAME
OUT_NAME = "ladrillo_table_a2.md"
DEFAULT_OUTDIR = os.path.join(REPO, "outputs", "tables")

ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
ap.add_argument("--pca", default=None, help=f"the PCA file (default <outdir>/{PCA_NAME})")
_a = ap.parse_args()
OUTDIR = os.path.abspath(_a.outdir)
PCA_FILE = os.path.abspath(_a.pca) if _a.pca else os.path.join(OUTDIR, PCA_NAME)

IDENT_MAX, PROP_MIN, LOAD_MIN = 0.20, 0.80, 0.30
NAME = {"ais_mu": "µ", "ais_bedheight0": "b₀", "ais_slope": "slope", "ais_iceflow0": "f₀", "ais_precip0_LOG": "ln P₀",
        "ais_runoff_Ton": "T_on", "ais_c": "c", "ais_gmst_amp": "amp", "antarctic_alpha": "α_DAIS", "antarctic_nu": "ν",
        "antarctic_temp_threshold": "T_crit", "anto_alpha": "a_ANTO", "anto_beta": "b_ANTO", "antarctic_lambda": "λ",
        "antarctic_gamma": "γ", "antarctic_kappa": "κ_DAIS", "ais_ocean_temperature₀": "T_oc,0"}
d = pd.read_csv(PCA_FILE)
vcols = [c for c in d.columns if c.startswith("v_")]
ident = d[d.ratio < IDENT_MAX].sort_values("ratio")
nprop = int((d.ratio > PROP_MIN).sum()); npart = len(d) - nprop - len(ident)
rows = []
for _, r in ident.iterrows():
    v = {c[2:]: r[c] for c in vcols}
    terms = sorted(((k, x) for k, x in v.items() if abs(x) >= LOAD_MIN), key=lambda t: -abs(t[1]))
    sgn = 1 if terms[0][1] > 0 else -1          # sign convention: leading loading positive
    combo = " ".join(f"{'+' if sgn*x > 0 else '−'} {abs(x):.2f} {NAME[k]}" for k, x in terms).lstrip("+ ")
    rows.append((int(r.pc), combo, 100 * r.ratio, r.rhat))
def lead(r, n=2):
    v = {c[2:]: r[c] for c in vcols}
    return ", ".join(NAME[k] for k, _ in sorted(v.items(), key=lambda t: -abs(t[1]))[:n])
prop = d[d.ratio > PROP_MIN].sort_values("pc")
prop_desc = "; ".join(f"{lead(r)}" for _, r in prop.iterrows())
cap = (f"**Table A2.** The directions in Antarctic parameter space that the observations identify. Principal components "
       f"of the {len(d)}-parameter Antarctic block of the posterior, each parameter standardised by its prior "
       f"standard deviation (the paleo-ensemble joint prior for the seven geometry parameters; Table A1's priors for the rest), "
       f"listed where the posterior variance along the component is below {int(IDENT_MAX*100)}% of the prior variance. "
       f"Loadings below {LOAD_MIN} are omitted; the sign is chosen so that the leading loading is positive. "
       f"Of the {len(d)} components, {nprop} {'retains' if nprop == 1 else 'retain'} more than {int(PROP_MIN*100)}% of {'its' if nprop == 1 else 'their'} prior variance (led by {prop_desc}) "
       f"and {npart} are partly identified; the 1900–2026 record therefore constrains {len(ident)} combinations of the "
       f"Antarctic parameters, not the parameters individually. R̂ is the split-R̂ of the component score across the four chains.")
lines = [cap, "", "| component | identified combination (loadings) | posterior variance, % of prior | R̂ |", "|---|---|---|---|"]
for pc, combo, pct, rh in rows:
    lines.append(f"| {pc} | {combo} | {pct:.1f} | {rh:.3f} |")
os.makedirs(OUTDIR, exist_ok=True)
out = os.path.join(OUTDIR, OUT_NAME)
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + f"\n\n<!-- Ladrillo {SCRIPT} | tag {TAG} | from {os.path.basename(PCA_FILE)} | "
            f"IDENT_MAX {IDENT_MAX} PROP_MIN {PROP_MIN} LOAD_MIN {LOAD_MIN} -->\n")
print("\n".join(lines)); print("wrote", out)
