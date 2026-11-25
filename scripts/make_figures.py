#!/usr/bin/env python
"""Generate the summary figures in docs/ from the cached data and benchmark JSON.

Figures (all from REAL computed numbers):
  1. docs/vs_enrichment_curves.png  - active-recovery curves, virtual screening
  2. docs/fp_metric_heatmap.png     - BEDROC across fingerprint x metric (VS)
  3. docs/task_winners.png          - best fingerprint+metric per benchmark task
"""
import json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)
import fsg
DATA, DOCS = os.path.join(ROOT, "data"), os.path.join(ROOT, "docs")
res = json.load(open(os.path.join(DOCS, "benchmark_results.json")))

act = fsg.load_mols(pd.read_csv(f"{DATA}/dude_egfr_actives.csv"))
dec = fsg.load_mols(pd.read_csv(f"{DATA}/dude_egfr_decoys.csv"))
FP_NAMES = fsg.FP_NAMES

# ---- Figure 1: active-recovery curves (single-query averaged) --------------- #
N_QUERY = 25
qrng = np.random.default_rng(7)
query_idx = sorted(qrng.choice(len(act), size=N_QUERY, replace=False).tolist())
grid = np.linspace(0, 1, 101)
show = ["ECFP4-2048", "FCFP4-2048", "AtomPair", "RDKit", "MACCS"]
plt.figure(figsize=(7, 6))
for name in show:
    afps = fsg.make_fps(name, list(act["mol"]))
    dfps = fsg.make_fps(name, list(dec["mol"]))
    curves = []
    for q in query_idx:
        cand_a = [afps[i] for i in range(len(afps)) if i != q]
        a_scores = fsg.bulk_similarity("tanimoto", afps[q], cand_a)
        d_scores = fsg.bulk_similarity("tanimoto", afps[q], dfps)
        scores = np.array(list(a_scores) + list(d_scores))
        labels = np.array([1] * len(a_scores) + [0] * len(d_scores))
        order = np.argsort(-scores, kind="mergesort")
        y = labels[order]
        recall = np.cumsum(y) / y.sum()
        x = np.arange(1, len(y) + 1) / len(y)
        curves.append(np.interp(grid, x, recall))
    mean_curve = np.mean(curves, axis=0)
    plt.plot(grid * 100, mean_curve * 100, lw=2, label=name)
plt.plot([0, 100], [0, 100], "k--", lw=1, label="random")
plt.xlabel("% of ranked library screened")
plt.ylabel("% of EGFR actives recovered")
plt.title("Virtual screening: active recovery (DUD-E EGFR)\nligand-based, Tanimoto, averaged over 25 query actives")
plt.legend(loc="lower right", fontsize=9)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{DOCS}/vs_enrichment_curves.png", dpi=130)
plt.close()
print("wrote vs_enrichment_curves.png")

# ---- Figure 2: BEDROC heatmap fingerprint x metric (VS) --------------------- #
metrics = ["tanimoto", "dice", "cosine", "tversky", "russel"]
b1 = res["benchmark1_virtual_screening"]["scores"]
M = np.full((len(FP_NAMES), len(metrics)), np.nan)
for i, fp in enumerate(FP_NAMES):
    for j, mt in enumerate(metrics):
        key = f"{fp}|{mt}"
        if key in b1:
            M[i, j] = b1[key]["BEDROC"]
plt.figure(figsize=(7, 7))
im = plt.imshow(M, cmap="viridis", aspect="auto", vmin=np.nanmin(M), vmax=np.nanmax(M))
plt.colorbar(im, label="BEDROC (alpha=20), higher = better")
plt.xticks(range(len(metrics)), [m.capitalize() for m in metrics])
plt.yticks(range(len(FP_NAMES)), FP_NAMES)
for i in range(len(FP_NAMES)):
    for j in range(len(metrics)):
        if not np.isnan(M[i, j]):
            plt.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                     color="white" if M[i, j] < (np.nanmin(M)+np.nanmax(M))/2 else "black", fontsize=8)
plt.title("Virtual screening BEDROC by fingerprint x metric\n(DUD-E EGFR; grey = metric undefined for count fp)")
plt.tight_layout()
plt.savefig(f"{DOCS}/fp_metric_heatmap.png", dpi=130)
plt.close()
print("wrote fp_metric_heatmap.png")

# ---- Figure 3: best combo per task ------------------------------------------ #
def best(d, key, rev=True):
    items = [(k, v[key]) for k, v in d.items() if v[key] == v[key]]
    return max(items, key=lambda t: t[1]) if rev else min(items, key=lambda t: t[1])

b2 = res["benchmark2_scaffold_hopping"]["scores"]
b4 = res["benchmark4_knn_qsar_esol"]["scores"]
tasks = ["VS\n(BEDROC)", "Scaffold hop\n(BEDROC)", "kNN QSAR\n(R2)"]
winners = [best(b1, "BEDROC"), best(b2, "BEDROC"), best(b4, "R2")]
vals = [w[1] for w in winners]
labels = [w[0] for w in winners]
plt.figure(figsize=(7, 5))
bars = plt.bar(tasks, vals, color=["#3b7dd8", "#2ca089", "#d8743b"])
for b, lab in zip(bars, labels):
    plt.text(b.get_x() + b.get_width()/2, b.get_height()/2, lab.replace("|", "\n+ "),
             ha="center", va="center", rotation=0, fontsize=9, color="white", fontweight="bold")
plt.ylabel("best score")
plt.title("Best fingerprint + metric per task (real benchmark winners)")
plt.tight_layout()
plt.savefig(f"{DOCS}/task_winners.png", dpi=130)
plt.close()
print("wrote task_winners.png")
