import scvelo as scv
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np

adata_split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317_bbknn.h5ad')
adata_split2_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317_bbknn.h5ad')
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317_bbknn.h5ad')

# compute cosine similarity
n_cell = adata_split1_seed317.layers["velocity"].shape[0]
adata_split1_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed317.layers['velocity'], nan=0)
adata_split2_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed317.layers['velocity'], nan=0)
cos_sim = cosine_similarity(adata_split1_seed317.layers["velocity_rmNA"],
                            adata_split2_seed317.layers["velocity_rmNA"])
adata_total.obs["cos_sim"] = [cos_sim[i,i] for i in range(0,n_cell)]
adata_total.obs["cos_sim"] = pd.DataFrame(adata_total.obs["cos_sim"])
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="cos_sim",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_seed317_cos_similarity.png")
print("seed 317 done!")

del adata_total.obs["cos_sim"]
adata_split1_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320_bbknn.h5ad')
adata_split2_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320_bbknn.h5ad')

# compute cosine similarity
n_cell = adata_split1_seed320.layers["velocity"].shape[0]
adata_split1_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed320.layers['velocity'], nan=0)
adata_split2_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed320.layers['velocity'], nan=0)
cos_sim = cosine_similarity(adata_split1_seed320.layers["velocity_rmNA"],
                            adata_split2_seed320.layers["velocity_rmNA"])
adata_total.obs["cos_sim"] = [cos_sim[i,i] for i in range(0,n_cell)]
adata_total.obs["cos_sim"] = pd.DataFrame(adata_total.obs["cos_sim"])
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="cos_sim",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_seed320_cos_similarity.png")
print("seed 320 done!")


