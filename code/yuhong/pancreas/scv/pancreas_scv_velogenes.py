import scvelo as scv
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# read split counts data
adata_split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed317_split1_seurat.h5ad')
adata_split2_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed317_split2_seurat.h5ad')

## process split1
scv.pp.normalize_per_cell(adata_split1_seed317)
scv.pp.log1p(adata_split1_seed317)
scv.pp.moments(adata_split1_seed317, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1_seed317, svd_solver="arpack")
sc.pp.neighbors(adata_split1_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_seed317)
scv.tl.recover_dynamics(adata_split1_seed317)
scv.tl.velocity(adata_split1_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split1_seed317)
## process split2
scv.pp.normalize_per_cell(adata_split2_seed317)
scv.pp.log1p(adata_split2_seed317)
scv.pp.moments(adata_split2_seed317, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2_seed317, svd_solver="arpack")
sc.pp.neighbors(adata_split2_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_seed317)
scv.tl.recover_dynamics(adata_split2_seed317)
scv.tl.velocity(adata_split2_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split2_seed317)

np.where(~np.isnan(adata_split1_seed317.layers['velocity'][0] + adata_split2_seed317.layers['velocity'][0]))

##################################
##################################
# read split counts data
adata_split1_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed320_split1_seurat.h5ad')
adata_split2_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed320_split2_seurat.h5ad')

## process split1
scv.pp.normalize_per_cell(adata_split1_seed320)
scv.pp.log1p(adata_split1_seed320)
scv.pp.moments(adata_split1_seed320, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1_seed320, svd_solver="arpack")
sc.pp.neighbors(adata_split1_seed320, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_seed320)
scv.tl.recover_dynamics(adata_split1_seed320)
scv.tl.velocity(adata_split1_seed320, mode="dynamical")
scv.tl.velocity_graph(adata_split1_seed320)

## process split2
scv.pp.normalize_per_cell(adata_split2_seed320)
scv.pp.log1p(adata_split2_seed320)
scv.pp.moments(adata_split2_seed320, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2_seed320, svd_solver="arpack")
sc.pp.neighbors(adata_split2_seed320, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_seed320)
scv.tl.recover_dynamics(adata_split2_seed320)
scv.tl.velocity(adata_split2_seed320, mode="dynamical")
scv.tl.velocity_graph(adata_split2_seed320)


np.where(~np.isnan(adata_split1_seed320.layers['velocity'][0] + adata_split2_seed320.layers['velocity'][0]))



