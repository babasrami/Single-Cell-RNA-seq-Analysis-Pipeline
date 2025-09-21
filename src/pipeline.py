#!/usr/bin/env python3
"""
End-to-end Scanpy pipeline for adrenal medulla scRNA-seq (GSE147821).

Implements:
1) Basic QC and filtering
2) Normalization and HVG selection
3) Cell-cycle scoring
4) PCA / neighbors / UMAP / Leiden clustering
5) Re-clustering of adrenal medulla populations (SCPs, Chromaffin, Sympathoblasts)
6) Simple pseudotime analysis with DPT + binned expression heatmaps

Usage
-----
# Create env and install deps
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Download raw data (2+ GB)
bash scripts/download_data.sh

# Run the pipeline (figures will render to screen)
python src/pipeline.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

RAW_DIR = "data/raw"

# Map GSM IDs to sample labels (as used in the notebook)
SAMPLE_MAP = {
    'GSM4446535': 'week8_001',
    'GSM4446536': 'week9_063',
    'GSM4446537': 'week6_088',
    'GSM4446538': 'week14_123',
    'GSM4446539': 'week12_124',
    'GSM4446540': 'week8_125',
    'GSM4446541': 'week9_005',
    'GSM4446542': 'week11_006',
    'GSM4446543': 'week9_007',
    'GSM4734601': 'week8_016',
    'GSM4734602': 'week9_031_paraganglia',
    'GSM4734603': 'week12_035',
    'GSM4734604': 'week12_036_extraadrenal'
}

S_PHASE_GENES = ['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2',
                 'MCM6','CDCA7','DTL','PRIM1','UHRF1','HELLS','RFC2','RPA2','NASP',
                 'RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2',
                 'RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
                 'USP1','CLSPN','POLA1']

G2M_GENES = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80',
             'CKS2','NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4',
             'CCNB2','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1',
             'KIF20B','HMMR','CDC20','TTK','CDC25C','KIF2C']

def read_raw_h5_matrices(raw_dir):
    files = [f for f in os.listdir(raw_dir) if f.endswith(".h5")]
    if not files:
        raise FileNotFoundError(f"No .h5 files found in {raw_dir}. Run scripts/download_data.sh first.")
    adatas = []
    for fp in files:
        path = os.path.join(raw_dir, fp)
        print(f"Loading: {path}")
        gsm_id = fp.split('_')[0]
        smp = SAMPLE_MAP.get(gsm_id, f"Sample_{gsm_id}")
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()
        adata.obs["sample_id"] = smp
        adatas.append(adata)
    return ad.concat(adatas)

def basic_qc(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # thresholds per the notebook
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    return adata

def normalize_and_hvg(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    return adata[:, adata.var.highly_variable]

def cell_cycle_score(adata):
    adata.var_names_make_unique()
    s_genes = [g for g in S_PHASE_GENES if g in adata.var_names]
    g2m_genes = [g for g in G2M_GENES if g in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata

def embed_and_cluster(adata):
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata

def recluster_adrenal(adata):
    # Choose cluster IDs of interest after inspecting UMAP (placeholder values)
    target_clusters = ['2','5','7']
    sub = adata[adata.obs['leiden'].isin(target_clusters)].copy()
    sc.pp.pca(sub); sc.pp.neighbors(sub); sc.tl.umap(sub); sc.tl.leiden(sub)
    # Map to labels (placeholder mapping—adjust after inspection)
    mapping = {'0':'SCPs','1':'Chromaffin cells','2':'Sympathoblasts'}
    sub.obs['cell_type_annotation'] = sub.obs['leiden'].map(mapping).astype('category')
    return sub

def simple_dpt_heatmap(adata, genes, n_bins=20, title="Trajectory"):
    # Neighborhood and pseudotime
    sc.pp.pca(adata); sc.pp.neighbors(adata)
    # diffmap improves DPT robustness
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    import seaborn as sns
    pt = adata.obs['dpt_pseudotime']
    order = np.argsort(pt)
    valid = [g for g in genes if g in adata.var_names]
    if not valid:
        print(f"No valid genes from {genes}; skipping {title}.")
        return
    df = adata.to_df()[valid].iloc[order]
    bins = pd.qcut(pt.iloc[order], q=n_bins, labels=False)
    df['bin'] = bins.values
    binned = df.groupby('bin').mean()
    ax = sns.heatmap(binned, cmap='viridis')
    ax.set_title(f"{title}: Gene Expression Dynamics")
    ax.set_xlabel("Genes"); ax.set_ylabel("Pseudotime Bins")
    plt.show()

def main():
    adata = read_raw_h5_matrices(RAW_DIR)
    print(adata)
    adata = basic_qc(adata)
    adata = normalize_and_hvg(adata)
    adata = cell_cycle_score(adata)
    adata = embed_and_cluster(adata)

    # Re-cluster adrenal medulla subset (optional)
    adrenal = recluster_adrenal(adata)
    sc.pl.umap(adrenal, color=['leiden','cell_type_annotation'], legend_loc='on data')

    # Simple trajectories (marker-focused heatmaps)
    scp = ['SOX10','PLP1','ISL1']
    chrom = ['CHGA','PNMT']
    symp = ['ELAVL4','PRPH']

    simple_dpt_heatmap(adata, list(set(scp+chrom)), title="SCP vs Chromaffin")
    simple_dpt_heatmap(adata, list(set(scp+symp)), title="SCP vs Sympathoblast")
    simple_dpt_heatmap(adata, list(set(chrom+symp)), title="Chromaffin vs Sympathoblast")

if __name__ == "__main__":
    main()
