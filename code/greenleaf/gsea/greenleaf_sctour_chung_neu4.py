# https://github.com/linnykos/veloUncertainty/blob/main/veloUncertainty/v4_functions_transMat.py

import scvelo as scv
import numpy as np
import collections

import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

import anndata as ad
import numpy as np
import scipy.sparse as sp
import mygene
import pandas as pd

import sys
sys.path.append('/home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/veloUncertainty')
from laplacian import *

# Load the AnnData object
adata_kernel = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/seed317/sct_GPC/adata_glf_sct_GPC_total_v4_outputAdded.h5ad")
adata_full = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/seed317/sct/adata_glf_sct_total_v4_outputAdded.h5ad")
adata_full, velocity_graph = prepare_adata_laplacian(adata_kernel=adata_kernel, adata_full=adata_full, subset_key='cluster_name', subset_value='excitatory neuron 4')
scores = compute_laplacian(adata=adata_full, velocity_graph=velocity_graph)

# np.nanquantile(scores, [0.0, 0.25, 0.5, 0.75, 1.0])
# np.isnan(scores).sum()

print("Finished computing all scores!", flush=True)

# Initialize
mg = mygene.MyGeneInfo()
# List of all 2000 gene IDs from your dataset
ensembl_ids = adata_full.var_names.tolist()  # Get all gene IDs from your AnnData object
# Query the mapping
out = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
# Build a mapping: Ensembl ID -> gene symbol
id_to_symbol = {entry['query']: entry.get('symbol', None) for entry in out}
# Print a few mappings to check
for k, v in list(id_to_symbol.items())[:10]:
    print(f"{k} → {v}")

# Count how many genes had no symbol
n_missing = sum(1 for symbol in id_to_symbol.values() if symbol is None)
# Total number of genes
n_total = len(id_to_symbol)
# Print summary
print(f"{n_missing} out of {n_total} Ensembl IDs did not have a gene symbol ({n_missing/n_total:.2%}).", flush=True)

# Get gene names from adata
gene_names = adata_full.var_names.tolist()
# Map Ensembl IDs to gene symbols
gene_symbols = [id_to_symbol.get(gene, None) for gene in gene_names]
# Create a DataFrame
scores_df = pd.DataFrame({
    'gene': gene_names,
    'symbol': gene_symbols,
    'score': scores
})

# Save to CSV
output_path = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/Writeup15_greenleaf_gene_laplacian_scores_sct_chung_neu4.csv"
scores_df.to_csv(output_path, index=False)

print(f"Saved scores to {output_path}", flush=True)