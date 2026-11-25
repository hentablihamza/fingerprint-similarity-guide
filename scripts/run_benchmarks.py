#!/usr/bin/env python
"""Run the four fingerprint x similarity benchmarks and save REAL results.

Writes docs/benchmark_results.json (consumed by Notebook 12 and the README) and
the summary figures in docs/. Every number is computed here from the cached data;
nothing is hard-coded. Deterministic (fixed seeds, full-data leave-one-out where
possible), so re-running reproduces the same tables.

Usage:  python scripts/run_benchmarks.py
"""
import json, os, sys, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)
import fsg
from rdkit.Chem import rdMolDescriptors
from rdkit import Chem

DATA = os.path.join(ROOT, "data")
DOCS = os.path.join(ROOT, "docs")
os.makedirs(DOCS, exist_ok=True)

def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# --------------------------------------------------------------------------- #
# Load data
# --------------------------------------------------------------------------- #
log("loading data ...")
act = fsg.load_mols(pd.read_csv(f"{DATA}/dude_egfr_actives.csv"))
dec = fsg.load_mols(pd.read_csv(f"{DATA}/dude_egfr_decoys.csv"))
egfr = fsg.load_mols(pd.read_csv(f"{DATA}/egfr_chembl37_ic50.csv"))
esol = fsg.load_mols(pd.read_csv(f"{DATA}/esol_delaney.csv"))
log(f"DUD-E actives={len(act)} decoys={len(dec)} | ChEMBL EGFR={len(egfr)} | ESOL={len(esol)}")

FP_NAMES = fsg.FP_NAMES
results = {"fingerprints": FP_NAMES}

# --------------------------------------------------------------------------- #
# Precompute fingerprints once per set
# --------------------------------------------------------------------------- #
def fps_for(name, mols):
    return fsg.make_fps(name, list(mols))

log("precomputing fingerprints (DUD-E) ...")
act_fps = {n: fps_for(n, act["mol"]) for n in FP_NAMES}
dec_fps = {n: fps_for(n, dec["mol"]) for n in FP_NAMES}

# --------------------------------------------------------------------------- #
# Benchmark 1: virtual screening recovery (ligand-based, single-query averaged)
#   Standard protocol: each of several seed actives is used ALONE as the query;
#   the remaining library (other actives + decoys) is ranked by similarity to it.
#   Metrics are averaged over the seeds. A single query (vs. max-over-all-542
#   analogues) avoids DUD-E analogue-bias saturation and so discriminates fps.
# --------------------------------------------------------------------------- #
log("Benchmark 1: virtual screening recovery ...")
N_QUERY = 25
qrng = np.random.default_rng(7)
query_idx = sorted(qrng.choice(len(act), size=N_QUERY, replace=False).tolist())

def avg_retrieval(name, metric, query_indices, scaffold_hop=False):
    return fsg.ligand_based_screen(
        act_fps[name], dec_fps[name], metric, query_indices,
        scaffolds=act_scaf if scaffold_hop else None, scaffold_hop=scaffold_hop)

# scaffolds are needed by both B1 (for the hop variant) and B2; compute now
from rdkit.Chem.Scaffolds import MurckoScaffold
def generic_scaffold(mol):
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        gen = MurckoScaffold.MakeScaffoldGeneric(core)
        return Chem.MolToSmiles(gen)
    except Exception:
        return None
act_scaf = [generic_scaffold(m) for m in act["mol"]]
n_scaf = len(set(s for s in act_scaf if s))

b1 = {}
for name in FP_NAMES:
    for metric in fsg.valid_metrics(name):
        b1[f"{name}|{metric}"] = avg_retrieval(name, metric, query_idx, scaffold_hop=False)
    log(f"  done {name}")
results["benchmark1_virtual_screening"] = {"n_queries": N_QUERY, "scores": b1}

# --------------------------------------------------------------------------- #
# Benchmark 2: scaffold hopping (single-query, different-scaffold recovery)
#   Same single-query screen as B1, but same-generic-scaffold actives are removed
#   so only recovery of actives on a DIFFERENT core counts. This is the scaffold
#   hop where exact ECFP tends to fall behind feature/reduced fingerprints.
# --------------------------------------------------------------------------- #
log("Benchmark 2: scaffold hopping ...")
log(f"  {n_scaf} distinct generic scaffolds among {len(act)} actives")
b2 = {}
for name in FP_NAMES:
    for metric in fsg.valid_metrics(name):
        b2[f"{name}|{metric}"] = avg_retrieval(name, metric, query_idx, scaffold_hop=True)
    log(f"  done {name}")

# Secondary, more sensitive scaffold-hop signal: how much active-active similarity
# survives when the two actives sit on DIFFERENT generic scaffolds. retention =
# mean(cross-scaffold sim) / mean(same-scaffold sim); higher = better at bridging
# scaffolds. Computed with Tanimoto so only the fingerprint varies.
log("Benchmark 2b: cross-scaffold similarity retention ...")
retention = {}
for name in FP_NAMES:
    afps = act_fps[name]
    same, cross = [], []
    for i in range(len(afps)):
        sims = fsg.bulk_similarity("tanimoto", afps[i], afps[i + 1:])
        for off, s in enumerate(sims):
            j = i + 1 + off
            (same if act_scaf[i] == act_scaf[j] else cross).append(s)
    ms = float(np.mean(same)) if same else float("nan")
    mc = float(np.mean(cross)) if cross else float("nan")
    retention[name] = {
        "same_scaffold_mean_sim": round(ms, 3),
        "cross_scaffold_mean_sim": round(mc, 3),
        "retention_ratio": round(mc / ms, 3) if ms else None,
    }
    log(f"  {name}: retention={retention[name]['retention_ratio']}")

results["benchmark2_scaffold_hopping"] = {
    "n_generic_scaffolds": n_scaf, "n_queries": N_QUERY,
    "scores": b2, "retention": retention,
}

# --------------------------------------------------------------------------- #
# Benchmark 3: activity-cliff sensitivity (ChEMBL EGFR pIC50, Tanimoto)
#   Among near-neighbour pairs (sim >= 0.7), what fraction are activity cliffs
#   (|dpIC50| >= 2)? And Spearman corr between similarity and -|dpIC50|.
# --------------------------------------------------------------------------- #
log("Benchmark 3: activity-cliff sensitivity ...")
rng = np.random.default_rng(42)
samp = egfr.sample(n=1500, random_state=42).reset_index(drop=True)
y = samp["p_activity"].to_numpy()
from scipy.stats import spearmanr
SIM_T, CLIFF_D = 0.70, 2.0
b3 = {}
for name in FP_NAMES:
    fps = fps_for(name, samp["mol"])
    sim_vals, dact_vals = [], []
    n_nn, n_cliff = 0, 0
    for i in range(len(fps)):
        sims = fsg.bulk_similarity("tanimoto", fps[i], fps[i + 1:])
        for off, s in enumerate(sims):
            j = i + 1 + off
            d = abs(y[i] - y[j])
            sim_vals.append(s); dact_vals.append(d)
            if s >= SIM_T:
                n_nn += 1
                if d >= CLIFF_D:
                    n_cliff += 1
    rho, _ = spearmanr(sim_vals, [-d for d in dact_vals])
    b3[name] = {
        "n_near_neighbor_pairs": int(n_nn),
        "n_cliff_pairs": int(n_cliff),
        "cliff_fraction": round(n_cliff / n_nn, 3) if n_nn else None,
        "spearman_sim_vs_activity_closeness": round(float(rho), 3),
    }
    log(f"  done {name}: nn={n_nn} cliffs={n_cliff}")
results["benchmark3_activity_cliffs"] = {
    "n_sampled": int(len(samp)), "sim_threshold": SIM_T, "cliff_delta": CLIFF_D,
    "metric": "tanimoto", "scores": b3,
}

# --------------------------------------------------------------------------- #
# Benchmark 4: nearest-neighbour QSAR regression (ESOL logS, 5-fold CV)
#   k-NN (k=5) with distance = 1 - similarity. Report R2 and RMSE.
# --------------------------------------------------------------------------- #
log("Benchmark 4: nearest-neighbour QSAR (ESOL) ...")
from sklearn.model_selection import KFold
yv = esol["logS"].to_numpy()
kf = KFold(n_splits=5, shuffle=True, random_state=42)
K = 5
b4 = {}
esol_fps_all = {n: fps_for(n, esol["mol"]) for n in FP_NAMES}
for name in FP_NAMES:
    fps = esol_fps_all[name]
    for metric in fsg.valid_metrics(name):
        preds = np.zeros(len(yv))
        for tr, te in kf.split(fps):
            tr_fps = [fps[i] for i in tr]
            for ti in te:
                sims = np.array(fsg.bulk_similarity(metric, fps[ti], tr_fps))
                top = np.argsort(-sims)[:K]
                w = sims[top]
                if w.sum() <= 0:
                    preds[ti] = yv[tr].mean()
                else:
                    preds[ti] = np.average(yv[tr][top], weights=w)
        ss_res = float(np.sum((yv - preds) ** 2))
        ss_tot = float(np.sum((yv - yv.mean()) ** 2))
        r2 = 1 - ss_res / ss_tot
        rmse = float(np.sqrt(np.mean((yv - preds) ** 2)))
        b4[f"{name}|{metric}"] = {"R2": round(r2, 3), "RMSE": round(rmse, 3)}
    log(f"  done {name}")
results["benchmark4_knn_qsar_esol"] = {"k": K, "cv_folds": 5, "endpoint": "logS", "scores": b4}

# --------------------------------------------------------------------------- #
# Save
# --------------------------------------------------------------------------- #
with open(os.path.join(DOCS, "benchmark_results.json"), "w") as f:
    json.dump(results, f, indent=2)
log("wrote docs/benchmark_results.json")
