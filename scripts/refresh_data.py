#!/usr/bin/env python
"""Refresh the cached datasets for fingerprint-similarity-guide.

Downloads the DUD-E EGFR virtual-screening benchmark (actives + property-matched
decoys) and subsamples the decoys to a size that runs quickly in Colab/CI. The
EGFR ChEMBL inhibitor set (with pIC50) and the Delaney ESOL solubility set are
reused from the author's other repositories and documented in data/README.md.

No measurements are invented. Decoy subsampling uses a fixed random seed so the
cached files are reproducible. Re-run with:  python scripts/refresh_data.py
"""
import json, os, random, sys, urllib.request
from datetime import date
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
DUDE = "https://dude.docking.org/targets/egfr"
N_DECOYS = 3000
SEED = 42

def canon(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return Chem.MolToSmiles(m)

def props(smiles):
    m = Chem.MolFromSmiles(smiles)
    return dict(
        mw=round(Descriptors.MolWt(m), 1),
        logp=round(Crippen.MolLogP(m), 2),
        hbd=rdMolDescriptors.CalcNumHBD(m),
        hba=rdMolDescriptors.CalcNumHBA(m),
        rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
    )

def fetch(name, local_hint=None):
    """Download a DUD-E file; fall back to a pre-downloaded copy in /tmp."""
    url = f"{DUDE}/{name}"
    try:
        with urllib.request.urlopen(url, timeout=180) as r:
            return r.read().decode("utf-8", "replace")
    except Exception as e:  # pragma: no cover - network
        if local_hint and os.path.exists(local_hint):
            print(f"  (network failed: {e}; using {local_hint})")
            return open(local_hint).read()
        raise

def main():
    os.makedirs(DATA, exist_ok=True)
    print("Fetching DUD-E EGFR actives ...")
    act_raw = fetch("actives_final.ism", "/tmp/dude_egfr_actives.ism")
    print("Fetching DUD-E EGFR decoys ...")
    dec_raw = fetch("decoys_final.ism", "/tmp/dude_egfr_decoys.ism")

    actives = []
    for line in act_raw.strip().splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        smi, zinc, chembl = parts[0], parts[1], parts[2]
        c = canon(smi)
        if c:
            actives.append((chembl, zinc, c))

    decoys_all = []
    for line in dec_raw.strip().splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        smi, zinc = parts[0], parts[1]
        decoys_all.append((zinc, smi))
    rng = random.Random(SEED)
    rng.shuffle(decoys_all)
    decoys = []
    for zinc, smi in decoys_all:
        c = canon(smi)
        if c:
            decoys.append((zinc, c))
        if len(decoys) >= N_DECOYS:
            break

    with open(os.path.join(DATA, "dude_egfr_actives.csv"), "w") as f:
        f.write("chembl_id,zinc_id,smiles,label\n")
        for chembl, zinc, smi in actives:
            f.write(f"{chembl},{zinc},{smi},1\n")
    with open(os.path.join(DATA, "dude_egfr_decoys.csv"), "w") as f:
        f.write("zinc_id,smiles,label\n")
        for zinc, smi in decoys:
            f.write(f"{zinc},{smi},0\n")

    # property summary to document decoy matching
    import statistics as st
    def summ(rows, idx):
        vals = {k: [] for k in ["mw", "logp", "hbd", "hba", "rotb"]}
        for r in rows:
            p = props(r[idx])
            for k in vals:
                vals[k].append(p[k])
        return {k: {"mean": round(st.mean(v), 2), "median": round(st.median(v), 2)}
                for k, v in vals.items()}

    prov = {
        "dude_egfr": {
            "retrieved": str(date.today()),
            "source": DUDE,
            "target": "EGFR (DUD-E target 'egfr', UniProt P00533)",
            "n_actives": len(actives),
            "n_decoys_total": len(decoys_all),
            "n_decoys_cached": len(decoys),
            "decoy_subsample_seed": SEED,
            "attribution": "DUD-E is provided free of charge by the Irwin and Shoichet Laboratories, UCSF.",
            "property_summary_actives": summ(actives, 2),
            "property_summary_decoys": summ(decoys, 1),
            "citation": "Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. "
                        "Directory of Useful Decoys, Enhanced (DUD-E). J. Med. Chem. 2012, 55, 6582-6594.",
        },
        "egfr_chembl37_ic50": {
            "source": "ChEMBL release 37, target CHEMBL203 (human EGFR), binding IC50",
            "curation": "author's chembl-easy pipeline",
            "n_compounds": 6871,
            "attribution": "ChEMBL data provided by EMBL-EBI under CC BY-SA 3.0",
        },
        "esol_delaney": {
            "source": "DeepChem/MoleculeNet distribution of the Delaney ESOL set",
            "n_compounds": 1128,
            "citation": "Delaney, J. S. J. Chem. Inf. Comput. Sci. 2004, 44, 1000-1005.",
        },
    }
    with open(os.path.join(DATA, "provenance.json"), "w") as f:
        json.dump(prov, f, indent=2)
    print(f"Wrote {len(actives)} actives, {len(decoys)} decoys, provenance.json")

if __name__ == "__main__":
    main()
