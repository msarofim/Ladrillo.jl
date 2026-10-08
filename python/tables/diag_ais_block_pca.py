#!/usr/bin/env python3
"""
diag_ais_block_pca.py -- principal components of the Antarctic block of the Ladrillo posterior, in PRIOR-sd
units, with the prior variance along each PC (identification) and a split R-hat across the four chains
(mixing). It is the input to Table A2 (ladrillo_table_a2.py): the unmixed direction is the PRIOR-DOMINATED one,
not an identified ridge.

  python python/tables/diag_ais_block_pca.py [--outdir=outputs/tables]

PORT of SLR-RFF-BRICK python/diag_ais_block_pca.py --tag=L27 --source=subsample, the basis of the paper's Table A2.
The numerical code is copied unchanged; only the I/O differs:
  posterior  data/posterior/ladrillo_L27_posterior.csv, the shipped 10k subsample (2nd halves of the four chains,
             every 400th iteration, in chain order: rows 1-2500 are chain 1, ...)
  priors     data/constants/paleo_geo_prior_ton.csv, the paleo joint prior of the seven geometry parameters, and
             data/constants/ladrillo_priors_L27.csv, the calibrator's own prior dump, for the other parameters
  output     <outdir>/diag_ais_block_pca.csv (+ provenance); the research file is
             outputs/diag_ais_block_pca_L27_sub10k.csv
The research script's --source=chains reads the MCMC chains, which this package does not ship; it is not ported.
python/tables/check_tables_identity.py gates the output against the research file.
"""
import argparse
import datetime
import os
import subprocess

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = "python/tables/diag_ais_block_pca.py"
TAG = "L27"                                       # the posterior's calibration
POST_FILE = os.path.join(REPO, "data", "posterior", f"ladrillo_{TAG}_posterior.csv")
GEO_PRIOR_FILE = os.path.join(REPO, "data", "constants", "paleo_geo_prior_ton.csv")
PRIOR_FILE = os.path.join(REPO, "data", "constants", f"ladrillo_priors_{TAG}.csv")
OUT_NAME = "diag_ais_block_pca.csv"
DEFAULT_OUTDIR = os.path.join(REPO, "outputs", "tables")

ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
OUTDIR = os.path.abspath(ap.parse_args().outdir)

SEEDS = [2026, 2027, 2028, 2029]
SUB_PER_CHAIN = 2500               # rows per chain, in chain order
SUB_FILE = POST_FILE
SRC_DESC = f"10k subsample {os.path.relpath(SUB_FILE, REPO)} ({len(SEEDS)} chains x {SUB_PER_CHAIN}, in chain order)"
TBAR = -17.992322907746065          # LADRILLO_TBAR_ANT (ladrillo_projection.jl); only used for a --precip-reparam chain
                                    # (== Ladrillo.tbar_ant() bit for bit, checked 2026-10-08)
GEO = ["ais_mu", "ais_bedheight0", "ais_slope", "ais_iceflow0", "ais_precip0_LOG", "ais_runoff_Ton", "ais_c"]
FD = ["ais_gmst_amp", "antarctic_alpha", "antarctic_nu", "antarctic_temp_threshold", "anto_alpha", "anto_beta",
      "antarctic_lambda", "antarctic_gamma", "antarctic_kappa", "ais_ocean_temperature₀"]
## tag-aware (2026-09-20, L27): parameters the tag does NOT sample (L27 holds lambda / T_crit / gamma at their
## paleo medians) are dropped from the block, read from the posterior's header -- never typed per tag.
_hdr0 = open(SUB_FILE, encoding="utf-8").readline().rstrip("\n").split(",")
FD = [c for c in FD if c in _hdr0]
USE = GEO + FD
NP = len(USE)
def _blocks():
    hdr = _hdr0
    cols = [c for c in USE if c in hdr] + (["ais_precip_u"] if "ais_precip_u" in hdr else [])
    d = pd.read_csv(SUB_FILE, usecols=cols)
    assert len(d) == len(SEEDS) * SUB_PER_CHAIN, f"{SUB_FILE}: {len(d)} rows, expected {len(SEEDS) * SUB_PER_CHAIN}"
    for k in range(len(SEEDS)):
        yield d.iloc[k * SUB_PER_CHAIN:(k + 1) * SUB_PER_CHAIN].copy()
chains = []
for d in _blocks():
    if "ais_precip_u" in d:
        d["ais_precip0_LOG"] = d.ais_precip_u - d.antarctic_kappa * TBAR
    chains.append(d[USE])
g = open(GEO_PRIOR_FILE, encoding="utf-8").read().splitlines()
rows = {l.split(",")[0]: l.split(",")[1:] for l in g if l and not l.startswith("#")}
sd = np.array([float(x) for x in rows["sd"]])
corr = np.array([[float(v) for v in l.split(",")[1:]] for l in g if l.startswith("corr")])
## the non-geometry priors are THIS TAG'S, read from the calibrator's own dump (--dump-priors)
pr = pd.read_csv(PRIOR_FILE).set_index("name")
psd = {c: float(pr.loc[c, "sigma"]) for c in FD}
pmu = {c: float(pr.loc[c, "mu"]) for c in FD}
S = np.zeros((NP, NP)); S[:7, :7] = np.outer(sd, sd) * corr
for i, c in enumerate(FD): S[7 + i, 7 + i] = psd[c] ** 2
mu = np.array([float(x) for x in rows["mean"]] + [pmu[c] for c in FD])
u = np.sqrt(np.diag(S))
Z = [(c.values - mu) / u for c in chains]
Cpost = np.cov(np.vstack(Z).T); Cpri = S / np.outer(u, u)
w, V = np.linalg.eigh(Cpost); o = np.argsort(w)[::-1]; w, V = w[o], V[:, o]
out = []
print(f"Antarctic block PCA, {TAG}, {SRC_DESC}, coordinates in prior-sd units; {len(Z)} chains x {len(Z[0])} draws")
print(f"{'PC':>3} {'post var':>9} {'prior var':>9} {'post/prior':>10} {'R-hat':>6}  loadings |>0.3|")
for k in range(NP):
    v = V[:, k]; pv = v @ Cpri @ v
    sc = [z @ v for z in Z]; m = np.array([x.mean() for x in sc]); vv = np.array([x.var(ddof=1) for x in sc]); n = len(sc[0])
    W = vv.mean(); B = n * m.var(ddof=1); rhat = np.sqrt(((n - 1) / n * W + B / n) / W)
    load = ", ".join(f"{USE[i]}{v[i]:+.2f}" for i in np.argsort(-abs(v)) if abs(v[i]) > 0.3)
    print(f"{k+1:>3} {w[k]:9.3f} {pv:9.3f} {w[k]/pv:10.2f} {rhat:6.3f}  {load}")
    out.append(dict(pc=k + 1, post_var=w[k], prior_var=pv, ratio=w[k] / pv, rhat=rhat, loadings=load,
                    **{f"v_{USE[i]}": v[i] for i in range(NP)}))
df = pd.DataFrame(out)
_commit = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
df["provenance"] = (f"Ladrillo {SCRIPT} | tag {TAG} | {SRC_DESC} | prior: paleo joint (geometry, "
                    f"{os.path.relpath(GEO_PRIOR_FILE, REPO)}), this tag's own priors for the rest "
                    f"({os.path.relpath(PRIOR_FILE, REPO)}) | split R-hat over 4 chains | commit {_commit} | {datetime.date.today()}")
os.makedirs(OUTDIR, exist_ok=True)
df.to_csv(os.path.join(OUTDIR, OUT_NAME), index=False)
print("cumulative variance first 3/5/8 PCs:", np.round(np.cumsum(w) / w.sum(), 2)[[2, 4, 7]])
print("wrote", os.path.join(OUTDIR, OUT_NAME))
