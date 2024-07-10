import datetime
import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/sct/" 

print_message_with_time("########### Start to read total, split1 and split2")
total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
adata_split1 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split1_allgenes.h5ad') # 9815 × 53801
adata_split2 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split2_allgenes.h5ad')

def scv_compute_velocity(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    ### batch correction
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("Batch correction done!")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

scv_compute_velocity(total) # running this
scv_compute_velocity(adata_split1)
scv_compute_velocity(adata_split2)

import numpy as np
common_genes_filter = np.intersect1d(np.array(adata_split1.var.index), np.array(adata_split2.var.index))
common_genes_filter.shape # 1499 common genes
np.intersect1d(np.array(adata_split1.var.index), np.array(total.var.index)).shape # 589
np.intersect1d(np.array(adata_split2.var.index), np.array(total.var.index)).shape # 580

#adata_split1.layers["velocity_rmNA"] = np.nan_to_num(adata_split1.layers['velocity'], nan=0)
#adata_split2.layers["velocity_rmNA"] = np.nan_to_num(adata_split2.layers['velocity'], nan=0)
np.sum(~np.isnan(adata_split1.layers['velocity'][0])) # Ngenes in split1 for velocity computation=9
np.sum(~np.isnan(adata_split2.layers['velocity'][0])) # Ngenes in split2 for velocity computation=9

velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
common_genes_velocity.shape # =1

# write data
total.write_h5ad(data_folder+'v2_erythroid/scv/adata_ery_scv_total_v2.h5ad')
adata_split1.write_h5ad(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split1_v2.h5ad')
adata_split2.write_h5ad(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split2_v2.h5ad')

import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
velo_df1 = pd.DataFrame(adata_split1.layers['velocity'], columns=adata_split1.var.index.tolist())
velo_df2 = pd.DataFrame(adata_split2.layers['velocity'], columns=adata_split2.var.index.tolist())
cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))

##### copied from an earlier script
scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)
scv.tl.score_genes_cell_cycle(adata)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)



