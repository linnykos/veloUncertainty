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

def chung_laplacian(A_csr):
    d_out = np.asarray(A_csr.sum(1)).ravel()
    P = sp.diags(1/np.maximum(d_out,1e-32)) @ A_csr      # row‑stochastic
    # power iteration for πᵀ
    pi = np.ones(A_csr.shape[0]) / A_csr.shape[0]
    for _ in range(10000):
        new = pi @ P
        if np.linalg.norm(new-pi,1) < 1e-12:
            break
        pi = new
    S  = sp.diags(np.sqrt(pi)) @ P @ sp.diags(1/np.sqrt(pi))
    H  = 0.5*(S + S.T)
    L  = sp.eye(A_csr.shape[0]) - H
    return L, pi


# Load the AnnData object
adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/seed317/velovi_woprep_GPC/adata_glf_velovi_woprep_GPC_total_GPC_outputAdded.h5ad")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

# Step 1: Check all values are non-negative
velocity_graph = vk.transition_matrix

adata2 = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/seed317/velovi_woprep/adata_glf_velovi_woprep_total_v4_outputAdded.h5ad")
# Get intersecting cell barcodes
shared_cells = np.intersect1d(adata.obs_names, adata2.obs_names)
# Subset both AnnData objects to only shared cells, and **sort them in the same order**
adata = adata[shared_cells].copy()
adata2 = adata2[shared_cells].copy()

# Subset
# Create the mask
keep_types = [
    "excitatory neuron 2",
    "excitatory neuron 4",
    "excitatory neuron 5",
    "excitatory neuron 6",
    "excitatory neuron 7"
]
mask = adata.obs['cluster_name'].isin(keep_types)
# Subset the AnnData object
adata = adata[mask].copy()
adata2 = adata2[mask].copy()
# Get indices of the cells you kept
indices = np.where(mask)[0]
# Subset velocity graph using those indices (preserve sparse structure)
velocity_graph = velocity_graph[indices, :][:, indices]

laplacian, pi = chung_laplacian(velocity_graph)

# Now loop over ALL genes (all columns of adata)
n_genes = adata2.n_vars  # This is 2000 genes in your case
scores = []

for gene_idx in range(n_genes):
    print(f"Working on gene {gene_idx + 1} out of {n_genes}", flush=True)
    
    # x : sparse column vector (CSR gives fast dot with L)
    x = adata2[:, gene_idx].X
    if not sp.issparse(x):
        x = sp.csr_matrix(x)          # handle the degenerate dense column
    else:
        x = x.tocsr()                 # make sure CSR for fast mat‑vec
    
    # ---- numerator  xᵀ L x  -------------------------------------------
    y      = laplacian @ x            # sparse (n×1)
    numer  = x.multiply(y).sum()      # element‑wise product, then sum
    
    # ---- π‑weighted mean / variance -----------------------------------
    idx_nz = x.indices                # positions of non‑zeros
    data   = x.data                   # the actual non‑zero values
    
    # mean_π = Σ π_i x_i  (zeros contribute nothing)
    mean_pi  = np.dot(pi[idx_nz], data)
    
    # variance: split into “non‑zero” and “zero” parts
    # (π sums to 1, so π_zero = 1 - Σₙₖ π_k where k are non‑zeros)
    pi_nz_sum   = pi[idx_nz].sum()
    pi_zero_sum = 1.0 - pi_nz_sum
    
    var_pi = np.dot(pi[idx_nz], (data - mean_pi) ** 2) + pi_zero_sum * (mean_pi ** 2)
    
    # guard against genes that are constant under π
    if var_pi < 1e-12:
        scores.append(np.nan)
        continue
    
    scores.append(float(numer / var_pi))

# Convert list to numpy array
scores = np.array(scores)

# np.nanquantile(scores, [0.0, 0.25, 0.5, 0.75, 1.0])
# np.isnan(scores).sum()

print("Finished computing all scores!", flush=True)

# Initialize
mg = mygene.MyGeneInfo()
# List of all 2000 gene IDs from your dataset
ensembl_ids = adata2.var_names.tolist()  # Get all gene IDs from your AnnData object
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
gene_names = adata2.var_names.tolist()
# Map Ensembl IDs to gene symbols
gene_symbols = [id_to_symbol.get(gene, None) for gene in gene_names]
# Create a DataFrame
scores_df = pd.DataFrame({
    'gene': gene_names,
    'symbol': gene_symbols,
    'score': scores
})

# Save to CSV
output_path = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/Writeup15_greenleaf_gene_laplacian_scores_velovi-woprep_chung_subset.csv"
scores_df.to_csv(output_path, index=False)

print(f"Saved scores to {output_path}", flush=True)