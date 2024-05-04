import scvelo as scv

adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/adata_pancreas_preprocess.h5ad")

## scatter plot without velocity arrows
scv.pl.scatter(adata, color='clusters', cmap='clusters_colors', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/pancreas_preprocess_scatter.png")

## velocity confidence
scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/pancreas_preprocess_velo_confidence.png")

## cosine similarity on preprocessed umap
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

# replace nan's to 0's in layers['velocity']
adata_split1_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed317.layers['velocity'], nan=0)
adata_split2_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed317.layers['velocity'], nan=0)
Ngenes_317s1 = np.sum(~np.isnan(adata_split1_seed317.layers['velocity'][0]))
Ngenes_317s2 = np.sum(~np.isnan(adata_split2_seed317.layers['velocity'][0]))

Ngenes_317common = np.sum(np.isnan(adata_split1_seed317.layers["velocity"][0] + adata_split2_seed317.layers["velocity"][0])==0)

# cosine similarity
cos_sim_seed317 = np.diag(cosine_similarity(adata_split1_seed317.layers["velocity_rmNA"], adata_split2_seed317.layers["velocity_rmNA"]))
adata.obs["cos_sim_seed317"] = cos_sim_seed317
adata.obs["cos_sim_seed317"] = pd.DataFrame(adata.obs["cos_sim_seed317"])
scv.pl.velocity_embedding_stream(adata, basis='umap',color="cos_sim_seed317",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/scvelo_seed317_cos_similarity_preumap.png")

## seed320
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
# replace nan's to 0's in layers['velocity']
adata_split1_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed320.layers['velocity'], nan=0)
adata_split2_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed320.layers['velocity'], nan=0)
Ngenes_320s1 = np.sum(~np.isnan(adata_split1_seed320.layers['velocity'][0]))
Ngenes_320s2 = np.sum(~np.isnan(adata_split2_seed320.layers['velocity'][0]))

Ngenes_320common = np.sum(np.isnan(adata_split1_seed320.layers["velocity"][0] + adata_split2_seed320.layers["velocity"][0])==0)

# cosine similarity
cos_sim_seed320 = np.diag(cosine_similarity(adata_split1_seed320.layers["velocity_rmNA"], adata_split2_seed320.layers["velocity_rmNA"]))

# add cosine similarities to total counts object
adata.obs["cos_sim_seed320"] = cos_sim_seed320
adata.obs["cos_sim_seed320"] = pd.DataFrame(adata.obs["cos_sim_seed320"])
scv.pl.velocity_embedding_stream(adata, basis='umap',color="cos_sim_seed320",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/scvelo_seed320_cos_similarity_preumap.png")
print("**************** seed320 cosine similarity plotted! ****************")






