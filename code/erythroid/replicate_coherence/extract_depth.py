import scanpy as sc
import scipy.io
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse

def compute_deseq2_dispersion_inputs(X, gene_names):
    """
    Prepare per-gene quantities for DESeq2-like dispersion shrinkage.
    
    Parameters
    ----------
    X : scipy.sparse._csr.csr_matrix
        Raw counts (cells x genes).
    alpha_hat : array-like, shape (n_genes,)
        Per-gene overdispersion estimates in DESeq2 parameterization
        (Var = mu + alpha * mu^2). Convert before calling if needed.
    design_rank : int, default=1
        Rank of the design matrix used when fitting dispersions.
        For an intercept-only model, design_rank = 1.
        Degrees of freedom will be n_cells - design_rank.
    
    Returns
    -------
    df_out : pandas.DataFrame
        Columns: "gene", "mean_norm", "alpha_hat".
    df : int
        Degrees of freedom used for dispersion estimation.
    """
    
    # Convert to dense numpy array
    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)
    
    # size factors (simple total-count)
    lib_size = X.sum(axis=1)
    median_lib = np.median(lib_size)
    if median_lib <= 0:
        raise ValueError("Median library size is non-positive; check counts.")
    
    size_factor = lib_size / median_lib
    
    # normalized counts and mean normalized expression per gene
    X_norm = X / size_factor[:, None]
    mean_norm = X_norm.mean(axis=0)
    
    df_out = pd.DataFrame({
        "gene": gene_names,
        "mean_norm": mean_norm
    })
    
    return df_out


data_folder = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

adata = ad.read_h5ad(data_folder+'Gastrulation/erythroid_lineage.h5ad')
spliced_matrix = adata.layers['spliced']
unspliced_matrix = adata.layers['unspliced']

###################


df_spliced = compute_deseq2_dispersion_inputs(X = spliced_matrix, )

df_out.to_csv("gene_dispersion_inputs.csv", index=False)

with open("df_for_dispersion.txt", "w") as f:
    f.write(str(df))
