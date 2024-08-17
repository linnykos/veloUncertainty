import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from sctour_misc import *

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
split_seed = 317
sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")


print_message_with_time("########### Read erythroid data")
adata = sc.read(data_folder+"Gastrulation/erythroid_lineage.h5ad")
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()
gene_names = adata.var.index.copy()

def train_sct_and_return_tnode(adata,n_top_genes,sct_seed=615):
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    adata.X = adata.X.astype(np.float32)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=n_top_genes, subset=True)
    tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
    return tnode

# preprocess
print_message_with_time("########### Training model for preprocessed adata")
tnode = train_sct_and_return_tnode(adata,n_top_genes=5000)
print_message_with_time("########### Computing velocities for preprocessed adata")
diff_mat_adata = compute_sctour_velocity(tnode, timestep=1/100)
adata.layers['velocity'] = diff_mat_adata

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
S_mat = S_mat[:,positions]
U_mat = U_mat[:,positions]
adata.layers['spliced_original'] = S_mat
adata.layers['unspliced_original'] = U_mat

print_message_with_time("########### Writing data for preprocessed adata")
adata.write(filename=data_folder+"v3_erythroid/sct/adata_ery_sct_preprocess_5kgenes.h5ad")
tnode.save_model(save_dir=data_folder+'v3_erythroid/sct/', save_prefix='tnode_ery_sct_preprocess_5kgenes')


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
    counts_adata1.write(data_folder+'v3_erythroid/sct/backup/counts_seed317_split1_5kgenes.h5ad')
    counts_adata2.write(data_folder+'v3_erythroid/sct/backup/counts_seed317_split2_5kgenes.h5ad')
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata_erythroid(split1[0],split1[1],total)
    adata2 = create_adata_erythroid(split2[0],split2[1],total)
    return adata1,adata2

adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=317)

adata_split1.write(data_folder+'v3_erythroid/sct/backup/seed317_split1_5kgenes.h5ad')
adata_split2.write(data_folder+'v3_erythroid/sct/backup/seed317_split2_5kgenes.h5ad')

# adata_split1 = create_adata_erythroid(S_split1,U_split1,adata)
print_message_with_time("########### Training model for split1")
tnode_split1 = train_sct_and_return_tnode(adata_split1,n_top_genes=2000)
print_message_with_time("########### Computing velocities for split1")
diff_mat_split1 = compute_sctour_velocity(tnode_split1, timestep=1/100)
adata_split1.layers['velocity'] = diff_mat_split1
print_message_with_time("########### Writing adata and tnode for split1")
adata_split1.write(data_folder+'v3_erythroid/sct/adata_ery_sct_seed317_split1.h5ad')
tnode_split1.save_model(save_dir=data_folder+'v3_erythroid/sct/', save_prefix='tnode_ery_sct_seed317_split1')

# adata_split2 = create_adata_erythroid(S_split2,U_split2,adata)
print_message_with_time("########### Training model for split2")
tnode_split2 = train_sct_and_return_tnode(adata_split2,n_top_genes=2000)
print_message_with_time("########### Computing velocities for split2")
diff_mat_split2 = compute_sctour_velocity(tnode_split2, timestep=1/100)
adata_split2.layers['velocity'] = diff_mat_split2
print_message_with_time("########### Writing adata and tnode for split2")
adata_split2.write(data_folder+'v3_erythroid/sct/adata_ery_sct_seed317_split2.h5ad')
tnode_split2.save_model(save_dir=data_folder+'v3_erythroid/sct/', save_prefix='tnode_ery_sct_seed317_split2')

print_message_with_time("########### Training model for total")
adata = sc.read(data_folder+"Gastrulation/erythroid_lineage.h5ad")
tnode_total = train_sct_and_return_tnode(adata,n_top_genes=2000)
print_message_with_time("########### Computing velocities for total")
diff_mat_total = compute_sctour_velocity(tnode_total, timestep=1/100)
adata.layers['velocity'] = diff_mat_total
print_message_with_time("########### Writing adata and tnode for total")
adata.write(data_folder+'v3_erythroid/sct/adata_ery_sct_total.h5ad')
tnode_total.save_model(save_dir=data_folder+'v3_erythroid/sct/', save_prefix='tnode_ery_sct_total')
print_message_with_time("########### All done")



