import datetime
import scvelo as scv
import scanpy as sc
import bbknn
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
split_seed = 317

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

print_message_with_time("########### Start to read total, split1 and split2")
adata = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()
gene_names = adata.var.index.copy()

print_message_with_time("########### Start to preprocess raw data, filter out 5000 genes, and get S_mat & U_mat")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=5000)

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
S_mat = S_mat[:,positions]
U_mat = U_mat[:,positions]
adata.layers['spliced_original'] = S_mat
adata.layers['unspliced_original'] = U_mat

def run_countsplit_with_overdispersion(S,U,split_seed):
    print_message_with_time("########### Estimating overdispersion parameters")
    overdisps_S = estimate_overdisps(S)
    overdisps_U = estimate_overdisps(U)
    print_message_with_time("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisps_S)
    u1, u2  = countsplit(U,overdisps=overdisps_U)
    return [[s1,u1],[s2,u2]]

def create_adata_erythroid(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.uns = {'celltype_colors':adata.uns['celltype_colors'].copy()}
    adata_split.obsm['X_pcaOriginal'] = adata_total.obsm['X_pca'].copy()
    adata_split.obsm['X_umapOriginal'] = adata_total.obsm['X_umap'].copy()
    adata_split.obs = pd.DataFrame(index=adata_total.obs.index)
    for obs_col in adata_total.obs.columns:
        adata_split.obs[obs_col] = adata_total.obs[obs_col].copy()
    adata_split.var = pd.DataFrame(index=adata_total.var.index)
    for var_col in adata_total.var.columns:
        adata_split.var[var_col] = adata_total.var[var_col].copy()
    return adata_split

def countsplit_and_create_adata(S,U,total,split_seed):
    print_message_with_time("########### Running the function for overdispersion estimation and countsplitting")
    split1,split2 = run_countsplit_with_overdispersion(S=S,U=U,split_seed=split_seed)
    print_message_with_time("########### Writing adata objects with counts only")
    counts_adata1 = ad.AnnData(X=split1[0].astype(np.float32))
    counts_adata1.layers["spliced"] = split1[0]
    counts_adata1.layers["unspliced"] = split1[1]
    counts_adata2 = ad.AnnData(X=split2[0].astype(np.float32))
    counts_adata2.layers["spliced"] = split2[0]
    counts_adata2.layers["unspliced"] = split2[1]
    counts_adata1.write(data_folder+'v3_erythroid/scv/backup/counts_seed317_split1_5kgenes.h5ad')
    counts_adata2.write(data_folder+'v3_erythroid/scv/backup/counts_seed317_split2_5kgenes.h5ad')
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata_erythroid(split1[0],split1[1],total)
    adata2 = create_adata_erythroid(split2[0],split2[1],total)
    return adata1,adata2

adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=317)

adata_split1.write(data_folder+'v3_erythroid/scv/backup/seed317_split1_5kgenes.h5ad')
adata_split2.write(data_folder+'v3_erythroid/scv/backup/seed317_split2_5kgenes.h5ad')

# adata.write(filename=data_folder+"v3_erythroid/scv/adata_ery_scv_preprocess_5kgenes.h5ad")

def scv_compute_velocity_erythroid(adata, n_top_genes=2000):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=n_top_genes)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    ### batch correction
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("********* Batch correction done! *********")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

print_message_with_time("########### Start to run scv procedure on splits and total")
scv_compute_velocity_erythroid(adata_split1)
scv_compute_velocity_erythroid(adata_split2)

total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
scv_compute_velocity_erythroid(total)

# write data
print_message_with_time("########### Start to write data")
total.write_h5ad(data_folder+'v3_erythroid/scv/adata_ery_scv_total_v3.h5ad')
adata_split1.write_h5ad(data_folder+'v3_erythroid/scv/adata_ery_scv_seed317_split1_v3.h5ad')
adata_split2.write_h5ad(data_folder+'v3_erythroid/scv/adata_ery_scv_seed317_split2_v3.h5ad')

print("************ Ngenes in being filtered out and N common genes ************")
common_genes_filter = np.intersect1d(np.array(adata_split1.var.index), np.array(adata_split2.var.index))
print(common_genes_filter.shape) # 1605 common genes
print(np.intersect1d(np.array(adata_split1.var.index), np.array(total.var.index)).shape) # 727
print(np.intersect1d(np.array(adata_split2.var.index), np.array(total.var.index)).shape) # 745

print("************ Ngenes in splits for velocity computation ************")
print(np.sum(~np.isnan(adata_split1.layers['velocity'][0]))) # Ngenes in split1 for velocity computation=9
print(np.sum(~np.isnan(adata_split2.layers['velocity'][0]))) # Ngenes in split2 for velocity computation=9

print("************ N common genes in splits for velocity computation ************")
velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
print(common_genes_velocity.shape) # 

print_message_with_time("########### All done")
