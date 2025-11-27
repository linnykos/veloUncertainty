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
    gene_names : array-like, shape (n_genes,)
        Gene names corresponding to X
    
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
gene_names = adata.var_names

spliced_matrix = adata.layers['spliced']
unspliced_matrix = adata.layers['unspliced']

###################

df_spliced = compute_deseq2_dispersion_inputs(X = spliced_matrix, gene_names = gene_names)
df_unspliced = compute_deseq2_dispersion_inputs(X = unspliced_matrix, gene_names = gene_names)

df_spliced.to_csv("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/erythroid_replicate_coherence/spliced_depth.csv", index=False)
df_unspliced.to_csv("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/erythroid_replicate_coherence/unspliced_depth.csv", index=False)
