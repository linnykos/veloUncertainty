import scvelo as scv
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np

adata_split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_split1_seurat.h5ad')
adata_split2_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_split2_seurat.h5ad')
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_total_seurat.h5ad')

scv.pp.normalize_per_cell(adata_split1_seed317)
scv.pp.log1p(adata_split1_seed317)
scv.pp.moments(adata_split1_seed317, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1_seed317, svd_solver="arpack")
sc.pp.neighbors(adata_split1_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_seed317)
scv.tl.recover_dynamics(adata_split1_seed317)
scv.tl.velocity(adata_split1_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split1_seed317)
scv.pl.velocity_embedding_stream(adata_split1_seed317, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo_seed317_split1_2000.png")
print("**************** seed317 split1 processed! ****************")


scv.pp.normalize_per_cell(adata_split2_seed317)
scv.pp.log1p(adata_split2_seed317)
scv.pp.moments(adata_split2_seed317, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2_seed317, svd_solver="arpack")
sc.pp.neighbors(adata_split2_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_seed317)
scv.tl.recover_dynamics(adata_split2_seed317)
scv.tl.velocity(adata_split2_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split2_seed317)
scv.pl.velocity_embedding_stream(adata_split2_seed317, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo_seed317_split2_2000.png")
print("**************** seed317 split2 processed! ****************")

# replace nan's to 0's in layers['velocity']
adata_split1_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed317.layers['velocity'], nan=0)
adata_split2_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed317.layers['velocity'], nan=0)
# cosine similarity
cos_sim = cosine_similarity(adata_split1_seed317.layers["velocity_rmNA"],
                            adata_split2_seed317.layers["velocity_rmNA"])
print("**************** seed317 cosine similarity computed! ****************")

# total counts data process
scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_total, svd_solver="arpack")
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.recover_dynamics(adata_total)
scv.tl.velocity(adata_total, mode="dynamical")
scv.tl.velocity_graph(adata_total)
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo_seed317_total_2000.png")
print("**************** total counts processed! ****************")

# add cosine similarities to total counts object
adata_total.obs["cos_sim_cell"] = np.diag(cos_sim)
adata_total.obs["cos_sim_cell"] = pd.DataFrame(adata_total.obs["cos_sim_cell"])
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="cos_sim_cell",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo_seed317_cos_similarity_2000.png")
print("**************** seed317 cosine similarity plotted! ****************")

# read split counts data
adata_split1_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/erythroid_seed320_split1_seurat.h5ad')
adata_split2_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/erythroid_seed320_split2_seurat.h5ad')
print("**************** read seed320 split counts! ****************")

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
scv.pl.velocity_embedding_stream(adata_split1_seed320, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo_seed320_split1_2000.png")
print("**************** seed320 split1 processed! ****************")

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
scv.pl.velocity_embedding_stream(adata_split2_seed320, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo_seed320_split2_2000.png")
print("**************** seed320 split2 processed! ****************")

# replace nan's to 0's in layers['velocity']
adata_split1_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed320.layers['velocity'], nan=0)
adata_split2_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed320.layers['velocity'], nan=0)
# cosine similarity
cos_sim = 0
cos_sim = cosine_similarity(adata_split1_seed320.layers["velocity_rmNA"],
                            adata_split2_seed320.layers["velocity_rmNA"])
print("**************** seed320 cosine similarity computed! ****************")

# add cosine similarities to total counts object
del adata_total.obs["cos_sim_cell"]
adata_total.obs["cos_sim_cell"] = np.diag(cos_sim)
adata_total.obs["cos_sim_cell"] = pd.DataFrame(adata_total.obs["cos_sim_cell"])
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="cos_sim_cell",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo_seed320_cos_similarity_2000.png")
print("**************** seed320 cosine similarity plotted! ****************")


