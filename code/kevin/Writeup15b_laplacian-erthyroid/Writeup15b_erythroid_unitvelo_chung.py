# https://github.com/linnykos/veloUncertainty/blob/main/veloUncertainty/v4_functions_transMat.py
# cellrank environment

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
adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_erythroid/seed317/utv_Mark/adata_ery_utv_Mark_total_outputAdded.h5ad")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

# Step 1: Check all values are non-negative
velocity_graph = vk.transition_matrix

adata2 = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_erythroid/seed317/utv/adata_ery_utv_total_v4_outputAdded.h5ad")
# Get intersecting cell barcodes
shared_cells = np.intersect1d(adata.obs_names, adata2.obs_names)
# Subset both AnnData objects to only shared cells, and **sort them in the same order**
adata = adata[shared_cells].copy()
adata2 = adata2[shared_cells].copy()

laplacian, pi = chung_laplacian(velocity_graph)

# Now loop over ALL genes (all columns of adata)
n_genes = adata2.n_vars  # This is 2000 genes in your case
scores = []

for gene_idx in range(n_genes):
    print(f"Working on gene {gene_idx + 1} out of {n_genes}", flush=True)
    
    # ---- numerator  xᵀ L x  -------------------------------------------
    x_dense = np.asarray(adata2[:, gene_idx].X).ravel()   # 1-D vector
    y       = laplacian.dot(x_dense)                      # L x
    numer   = x_dense.dot(y)                              # xᵀ L x

    nonzero_mask = x_dense != 0

    # mean_π = Σ π_i x_i  (zeros contribute nothing)
    mean_pi      = np.dot(pi[nonzero_mask], x_dense[nonzero_mask])

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

# Get gene names from adata
gene_names = adata2.var_names.tolist()
# Create a DataFrame
scores_df = pd.DataFrame({
    'gene': gene_names,
    'score': scores
})

# Save to CSV
output_path = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15b/Writeup15b_erythroid_gene_laplacian_scores_utv_chung.csv"
scores_df.to_csv(output_path, index=False)

print(f"Saved scores to {output_path}", flush=True)