# Datasets

All datasets here are small, public, and committed so every notebook runs offline. No measurements
were invented, averaged, or altered. Machine-readable provenance is in [`provenance.json`](provenance.json),
regenerated together with the cached files by [`../scripts/refresh_data.py`](../scripts/refresh_data.py).

## `dude_egfr_actives.csv` / `dude_egfr_decoys.csv` — virtual-screening benchmark

Used in notebooks 04 (virtual screening) and 05 (scaffold hopping).

- **542 actives** (`chembl_id`, `zinc_id`, `smiles`, `label=1`) and **3,000 property-matched decoys**
  (`zinc_id`, `smiles`, `label=0`), subsampled with a fixed seed (42) from the full 35,050 DUD-E
  decoys so the notebooks run quickly.
- DUD-E decoys are chosen to match the actives on molecular weight, logP, H-bond donors/acceptors,
  rotatable bonds, and net charge, while being topologically dissimilar — i.e. presumed inactive.
  The property summaries in `provenance.json` document the match.

**Source:** DUD-E, target `egfr` (UniProt P00533), https://dude.docking.org/targets/egfr. DUD-E is
provided free of charge by the Irwin and Shoichet Laboratories, UCSF.

**Cite:** Mysinger, M. M.; Carchia, M.; Irwin, J. J.; Shoichet, B. K. Directory of Useful Decoys,
Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking. *J. Med. Chem.* **2012**,
*55*, 6582–6594. https://doi.org/10.1021/jm300687e

## `egfr_chembl37_ic50.csv` — EGFR inhibitors with pIC50 (6,871 compounds)

Used in notebook 06 (activity cliffs) and 09/10 (substructure, clustering).

- `compound_id` — ChEMBL molecule ID
- `smiles` — curated structure
- `p_activity` — aggregated pIC50 (−log10 IC50 in mol/L; 6 = 1 µM, 9 = 1 nM)
- `active` — 1 if pIC50 ≥ 7, else 0

**Source:** ChEMBL release 37 (target CHEMBL203, human EGFR, binding IC50), curated with the author's
[chembl-easy](https://github.com/hentablihamza/chembl-easy) pipeline. Reused from the author's
[practical-cheminformatics-course](https://github.com/hentablihamza/practical-cheminformatics-course).

## `esol_delaney.csv` — aqueous solubility (1,128 compounds)

Used in notebook 07 (nearest-neighbour QSAR).

- `compound_id` — source compound name
- `smiles` — structure
- `logS` — measured aqueous solubility, log10(mol/L)

**Source:** the DeepChem / MoleculeNet distribution of the Delaney ESOL set, reused via the author's
[qsar-split-bench](https://github.com/hentablihamza/qsar-split-bench) repository.

**Cite:** Delaney, J. S. ESOL: Estimating Aqueous Solubility Directly from Molecular Structure.
*J. Chem. Inf. Comput. Sci.* **2004**, *44*, 1000–1005. https://doi.org/10.1021/ci034243x

## ChEMBL data attribution

ChEMBL data are provided by EMBL-EBI under the Creative Commons Attribution-ShareAlike 3.0 licence
(CC BY-SA 3.0). The ChEMBL-derived files keep their ChEMBL IDs. If you use them, please cite:

> Zdrazil, B. *et al.* The ChEMBL Database in 2023: a drug discovery platform spanning multiple
> bioactivity data types and time periods. *Nucleic Acids Res.* **2024**, *52*, D1180–D1192.
> https://doi.org/10.1093/nar/gkad1004
