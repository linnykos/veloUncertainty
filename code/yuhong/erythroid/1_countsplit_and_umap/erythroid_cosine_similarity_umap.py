import scvelo as scv
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
from matplotlib import colormaps

# process split counts data
adata_split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317.h5ad')
adata_split2_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317.h5ad')

scv.pp.normalize_per_cell(adata_split1_seed317)
scv.pp.log1p(adata_split1_seed317)
scv.pp.moments(adata_split1_seed317, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1_seed317, svd_solver="arpack")
sc.pp.neighbors(adata_split1_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_seed317)
scv.tl.recover_dynamics(adata_split1_seed317)
scv.tl.velocity(adata_split1_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split1_seed317)

scv.pp.normalize_per_cell(adata_split2_seed317)
scv.pp.log1p(adata_split2_seed317)
scv.pp.moments(adata_split2_seed317, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2_seed317, svd_solver="arpack")
sc.pp.neighbors(adata_split2_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_seed317)
scv.tl.recover_dynamics(adata_split2_seed317)
scv.tl.velocity(adata_split2_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split2_seed317)

# process total counts data
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317.h5ad')
scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_total, svd_solver="arpack")
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.recover_dynamics(adata_total)
scv.tl.velocity(adata_total, mode="dynamical")
scv.tl.velocity_graph(adata_total)

# compute cosine similarity
cos_sim = cosine_similarity(adata_split1_seed317.layers["velocity"],adata_split2_seed317.layers["velocity"])
adata_total.obs["cos_sim_cell"] = [cos_sim[i,i] for i in range(0,9815)]
adata_total.obs["cos_sim_cell"] = pd.DataFrame(adata_total.obs["cos_sim_cell"])
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="cos_sim_cell",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_cosine_similarity_seed317.png")



