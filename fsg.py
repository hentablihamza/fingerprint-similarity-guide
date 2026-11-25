"""fsg — shared helpers for the fingerprint-similarity-guide course.

A single, tested place for:
  * a **fingerprint factory** (Morgan/ECFP, FCFP, RDKit, atom pairs, torsions,
    MACCS, Avalon — binary and count flavours), and
  * **similarity metrics** (Tanimoto, Dice, Cosine, Tversky, Russel-Rao), and
  * small **retrieval-evaluation** functions (enrichment factor, BEDROC, ROC-AUC).

Every notebook imports these so the numbers are computed one way, everywhere.
Nothing here invents data: all functions operate on molecules you pass in.
"""
from __future__ import annotations
import math
import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator
from rdkit.Avalon import pyAvalonTools

RDLogger.DisableLog("rdApp.*")

# --------------------------------------------------------------------------- #
# Fingerprint factory
# --------------------------------------------------------------------------- #
# Generators are built once and reused (the RDKit-recommended pattern).
_morgan2_1024 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
_morgan2_2048 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_morgan3_2048 = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
_rdkitfp = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
_apgen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
_ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
# FCFP = feature-based Morgan (pharmacophore-flavoured atom invariants)
_feat_inv = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
_fcfp2_2048 = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048, atomInvariantsGenerator=_feat_inv)

# name -> (callable(mol) -> fingerprint, is_count, kind)
_FP_BUILDERS = {
    "ECFP4-1024":  (lambda m: _morgan2_1024.GetFingerprint(m), False, "binary"),
    "ECFP4-2048":  (lambda m: _morgan2_2048.GetFingerprint(m), False, "binary"),
    "ECFP6-2048":  (lambda m: _morgan3_2048.GetFingerprint(m), False, "binary"),
    "ECFP4-counts":(lambda m: _morgan2_2048.GetCountFingerprint(m), True, "count"),
    "FCFP4-2048":  (lambda m: _fcfp2_2048.GetFingerprint(m), False, "binary"),
    "RDKit":       (lambda m: _rdkitfp.GetFingerprint(m), False, "binary"),
    "AtomPair":    (lambda m: _apgen.GetFingerprint(m), False, "binary"),
    "Torsion":     (lambda m: _ttgen.GetFingerprint(m), False, "binary"),
    "MACCS":       (lambda m: MACCSkeys.GenMACCSKeys(m), False, "binary"),
    "Avalon":      (lambda m: pyAvalonTools.GetAvalonFP(m, nBits=2048), False, "binary"),
}

#: The default comparison set used across the course.
FP_NAMES = list(_FP_BUILDERS.keys())

def fp_is_count(name: str) -> bool:
    return _FP_BUILDERS[name][1]

def make_fp(name: str, mol):
    """One fingerprint for one molecule."""
    return _FP_BUILDERS[name][0](mol)

def valid_metrics(name: str):
    """Metrics that are defined for this fingerprint's representation.

    Count fingerprints use only the count-aware metrics (Tanimoto, Dice); the
    remaining set-overlap metrics are defined on plain bit vectors.
    """
    if fp_is_count(name):
        return ["tanimoto", "dice"]
    return ["tanimoto", "dice", "cosine", "tversky", "russel"]

def make_fps(name: str, mols):
    """A list of fingerprints for a list of RDKit molecules (skips None)."""
    build = _FP_BUILDERS[name][0]
    return [build(m) for m in mols if m is not None]

# --------------------------------------------------------------------------- #
# Similarity metrics
# --------------------------------------------------------------------------- #
# Tversky is asymmetric: alpha weights the query's bits, beta the target's.
# alpha=beta=1 -> Tanimoto; alpha=beta=0.5 -> Dice. For substructure-style
# "is the query contained in the target?" use alpha high, beta low (e.g. 0.9/0.1).
def bulk_similarity(metric: str, query_fp, target_fps, alpha=0.9, beta=0.1):
    """Similarity of one query fingerprint against many targets -> list[float].

    metric in {tanimoto, dice, cosine, tversky, russel}. Russel-Rao needs plain
    bit vectors (it counts shared on-bits / total bits), so it is not defined for
    count fingerprints.
    """
    m = metric.lower()
    if m == "tanimoto":
        return DataStructs.BulkTanimotoSimilarity(query_fp, target_fps)
    if m == "dice":
        return DataStructs.BulkDiceSimilarity(query_fp, target_fps)
    if m == "cosine":
        return DataStructs.BulkCosineSimilarity(query_fp, target_fps)
    if m == "russel":
        return DataStructs.BulkRusselSimilarity(query_fp, target_fps)
    if m == "tversky":
        return DataStructs.BulkTverskySimilarity(query_fp, target_fps, alpha, beta)
    raise ValueError(f"unknown metric: {metric}")

def similarity(metric: str, fp1, fp2, alpha=0.9, beta=0.1):
    """Similarity between two fingerprints (scalar)."""
    return bulk_similarity(metric, fp1, [fp2], alpha, beta)[0]

# --------------------------------------------------------------------------- #
# Retrieval evaluation (ranked list of actives vs decoys)
# --------------------------------------------------------------------------- #
# These wrap RDKit's authoritative rdkit.ML.Scoring.Scoring implementation so the
# enrichment, BEDROC, and AUC numbers match the published definitions exactly.
from rdkit.ML.Scoring import Scoring as _Scoring

def _ranked(scores, labels):
    """Pack (score, label) pairs sorted by score descending for RDKit Scoring."""
    pairs = sorted(zip(scores, labels), key=lambda t: -t[0])
    return [[p[1]] for p in pairs]  # col 0 holds the label

def roc_auc(scores, labels):
    """Area under the ROC curve (higher score = more active-like)."""
    return float(_Scoring.CalcAUC(_ranked(scores, labels), 0))

def enrichment_factor(scores, labels, fraction):
    """Enrichment factor at the top `fraction` (e.g. 0.01 for EF1%)."""
    return float(_Scoring.CalcEnrichment(_ranked(scores, labels), 0, [fraction])[0])

def bedroc(scores, labels, alpha=20.0):
    """BEDROC (Truchon & Bayly 2007); alpha=20 weights roughly the top ~8%."""
    return float(_Scoring.CalcBEDROC(_ranked(scores, labels), 0, alpha))

def screen_metrics(scores, labels):
    """EF1%, EF5%, BEDROC(20) and AUC for one ranked list -> dict."""
    return {
        "EF1%": round(enrichment_factor(scores, labels, 0.01), 2),
        "EF5%": round(enrichment_factor(scores, labels, 0.05), 2),
        "BEDROC": round(bedroc(scores, labels, 20.0), 3),
        "AUC": round(roc_auc(scores, labels), 3),
    }

def ligand_based_screen(active_fps, decoy_fps, metric, query_indices,
                        scaffolds=None, scaffold_hop=False):
    """Average single-query ligand-based screen used by Benchmarks 1 and 2.

    Each query index q uses active_fps[q] alone as the query; the remaining
    actives (optionally only those on a DIFFERENT generic scaffold) are the hits
    to recover and all decoys are the background. Returns metrics averaged over
    the query set. Identical code path for the notebooks and scripts/run_benchmarks.py
    so the numbers match exactly.
    """
    import numpy as _np
    agg = {"EF1%": [], "EF5%": [], "BEDROC": [], "AUC": []}
    for q in query_indices:
        qfp = active_fps[q]
        cand = []
        for i in range(len(active_fps)):
            if i == q:
                continue
            if scaffold_hop and scaffolds is not None and scaffolds[i] == scaffolds[q]:
                continue
            cand.append(active_fps[i])
        if not cand:
            continue
        a_scores = bulk_similarity(metric, qfp, cand)
        d_scores = bulk_similarity(metric, qfp, decoy_fps)
        scores = list(a_scores) + list(d_scores)
        labels = [1] * len(a_scores) + [0] * len(d_scores)
        m = screen_metrics(scores, labels)
        for k in agg:
            agg[k].append(m[k])
    return {k: round(float(_np.mean(v)), 3) for k, v in agg.items()}

def load_mols(df, smiles_col="smiles"):
    """Add a 'mol' column; drop rows RDKit cannot parse. Returns a copy."""
    out = df.copy()
    out["mol"] = out[smiles_col].apply(Chem.MolFromSmiles)
    return out[out["mol"].notnull()].reset_index(drop=True)
