import scvelo as scv
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bbknn

adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_total_seurat.h5ad')

scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
## batch correction
bbknn.bbknn(adata_total, batch_key='sequencing.batch')
adata_total.X = adata_total.X.toarray()
bbknn.ridge_regression(adata_total, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_total)
bbknn.bbknn(adata_total, batch_key='sequencing.batch')

sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.recover_dynamics(adata_total)
scv.tl.velocity(adata_total, mode="dynamical")
scv.tl.velocity_graph(adata_total)

scv.tl.velocity_confidence(adata_total)
scv.pl.scatter(adata_total, c='velocity_confidence', cmap='coolwarm', perc=[5, 95],
               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_total_confidence.png")





