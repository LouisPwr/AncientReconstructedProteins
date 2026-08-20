# Ancient Reconstructed Proteins — Supplementary Scripts

Scripts and notebooks supporting the manuscript **"Ancient Reconstructed Proteins: A Framework for Resurrecting Protein Structures from Million-Year-Old Metagenomes"**.

Contains the assembly/annotation Snakemake pipeline and the downstream plotting and analysis notebooks used to generate the results and figures in the paper.

## Contents

- `snakemake_pipeline/` — assembly (MEGAHIT, CarpeDeam), read mapping, MMseqs2, and Prokka annotation pipeline
- `01_ancientorfs_paths.ipynb` – `07_annotations.ipynb` — sequential downstream analysis notebooks (target selection, pathway matching, Foldseek structural search, annotation)
- `annotations.json`, `sample_id_greenland.txt` — supporting metadata

## Data

Sequence, structure, and annotation data referenced by these scripts are archived on Zenodo: [10.5281/zenodo.22027760](https://doi.org/10.5281/zenodo.22027760).
