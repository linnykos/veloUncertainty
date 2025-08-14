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

import sys
sys.path.append('/home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/veloUncertainty')
from laplacian import *

# Load the AnnData object
adata_kernel = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_pancreas/seed317/scv_Mark227/adata_pan_scv_Mark227_total.h5ad")
adata_full = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_pancreas/seed317/scv/adata_pan_scv_total_v4.h5ad")
adata_full, velocity_graph = prepare_adata_laplacian(adata_kernel=adata_kernel, adata_full=adata_full)
scores = compute_laplacian(adata=adata_full, velocity_graph=velocity_graph)

# np.nanquantile(scores, [0.0, 0.25, 0.5, 0.75, 1.0])
# np.isnan(scores).sum()

print("Finished computing all scores!", flush=True)

# Get gene names from adata
gene_names = adata_full.var_names.tolist()
# Create a DataFrame
scores_df = pd.DataFrame({
    'gene': gene_names,
    'score': scores
})


# Save to CSV
output_path = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15c/Writeup15c_pancreas_gene_laplacian_scores_scv_chung.csv"
scores_df.to_csv(output_path, index=False)

print(f"Saved scores to {output_path}", flush=True)