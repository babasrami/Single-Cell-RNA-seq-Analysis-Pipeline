# Adrenal scRNA Trajectory Analysis — GSE147821

End-to-end pipeline for single-cell RNA-seq (scRNA-seq) analysis using Scanpy. It takes raw 10x Genomics HDF5 files and walks you through cleaning the data (quality control, normalization), finding informative genes, reducing dimensions, clustering cells, visualizing them with UMAP, and arranging cells along a simple pseudotime trajectory to explore lineage or state changes. It’s meant as a clear starting template for students and researchers new to single-cell analysis. The repository includes a Jupyter notebook, a Python script, and pinned environments. You can run the included example dataset (GSE147821, human adrenal medulla) or swap in your own 10x files. Outputs include QC summaries, UMAP cluster plots, and basic gene-trend heatmaps across pseudotime.

## Contents
- `notebooks/` — Original Jupyter notebook from the task
- `src/pipeline.py` — End-to-end Scanpy pipeline (QC → HVG → clustering → re-clustering → DPT heatmaps)
- `scripts/download_data.sh` — GEO data downloader (~2 GB)
- `docs/Report.pdf` — Methods + interpretation
- `results/Code_Plots.pdf` — Saved outputs
- `requirements.txt` (pinned), `environment.yml` (conda), `LICENSE`, `.gitignore`

## Reproducible setup
**Conda (recommended)**
```bash
conda env create -f environment.yml
conda activate adrenal-scrna
```

**pip / venv**
```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Data
- GEO: **GSE147821** — human fetal adrenal medulla 10x matrices (multiple GSMs).
- Download raw matrices:
```bash
bash scripts/download_data.sh
```

## Run
```bash
python src/pipeline.py
```
Figures render to screen; to save, add `plt.savefig("results/fig.png", dpi=300, bbox_inches="tight")` before `plt.show()`.

## Notes
- Memory: ≥16–32 GB RAM recommended for local runs.
- The cluster→cell type mapping in `recluster_adrenal()` is a placeholder; inspect UMAPs and adjust for your data slice.

## Citations
- **Scanpy** — Wolf FA, Angerer P, Theis FJ. *Genome Biology* (2018) 19:15.
- **UMAP** — McInnes L, Healy J, Melville J. arXiv:1802.03426.
- **Leiden** — Traag VA, Waltman L, van Eck NJ. *Sci. Rep.* (2019) 9:5233.
- **Anndata** — Virshup I et al. *Genome Biology* (2021) 22:286.
- **Dataset** — NCBI GEO accession **GSE147821** (adrenal medulla scRNA-seq).

## License
MIT — see `LICENSE`.
