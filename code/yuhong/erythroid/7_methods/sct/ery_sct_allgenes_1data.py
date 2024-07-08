import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
#import scvelo as scv
import datetime


import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from sctour_misc import *

sct_seed = 615
split_seed = 317

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

# Example usage
print_message_with_time("Starting countsplit")
# time: start of script, countsplit, overdispersion, velocity method, cosine similarity, entire script

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/sct/" 


adata = sc.read(data_folder+"Gastrulation/erythroid_lineage.h5ad")
# countsplit
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("########### Estimating overdispersion paramter for spliced matrix ###########")
overdisp_S = estimate_overdisps(S_mat)
print_message_with_time("########### Estimating overdispersion paramter for unspliced matrix ###########")
overdisp_U = estimate_overdisps(U_mat)

np.random.seed(split_seed)
print_message_with_time("########### Countsplitting for spliced matrix ###########")
S_split1, S_split2  = countsplit(S_mat,overdisps=overdisp_S)
print_message_with_time("########### Countsplitting for unspliced matrix ###########")
U_split1, U_split2  = countsplit(U_mat,overdisps=overdisp_U)

# compute velocity
def create_adata_erythroid(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split)
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.obs['celltype'] = adata_total.obs['celltype'].copy()
    adata_split.obs['sequencing.batch'] = adata_total.obs['sequencing.batch'].copy()
    adata_split.var['Accession'] = adata_total.var['Accession'].index.copy()
    adata_split.uns['celltype_colors'] = adata_total.uns['celltype_colors'].copy()
    adata_split.obsm['X_pcaOriginal'] = adata_total.obsm['X_pca'].copy()
    adata_split.obsm['X_umapOriginal'] = adata_total.obsm['X_umap'].copy()
    return adata_split

def train_sct_and_return_tnode(adata):
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    adata.X = adata.X.astype(np.float32)
    #adata.layers['spliced'] = adata.layers['spliced'].astype(np.float32)
    #adata.layers['unspliced'] = adata.layers['unspliced'].astype(np.float32)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)
    tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
    return tnode

adata_split1 = create_adata_erythroid(S_split1,U_split1,adata)
adata_split1.write(data_folder+'v2_erythroid/adata_allgenes_seed317_split1.h5ad')
print_message_with_time("########### Split1 adata object created, start training model ###########")
tnode_split1 = train_sct_and_return_tnode(adata_split1)
print_message_with_time("########### Split1 model training completed, strat computing velocity ###########")
diff_mat_split1 = compute_sctour_velocity(tnode_split1, timestep=1/100)
print_message_with_time("########### Split1 velocity computed, start writing data ###########")
adata_split1.layers['velocity'] = diff_mat_split1
adata_split1.write(data_folder+'v2_erythroid/sct/adata_seed317_split1.h5ad')
tnode_split1.save_model(save_dir=data_folder+'v2_erythroid/sct/', save_prefix='tnode_seed317_split1')
print_message_with_time("########### Split1 data wrote ###########")

adata_split2 = create_adata_erythroid(S_split2,U_split2,adata)
adata_split2.write(data_folder+'v2_erythroid/adata_allgenes_seed317_split2.h5ad')
print_message_with_time("########### Split2 adata object created, start training model ###########")
tnode_split2 = train_sct_and_return_tnode(adata_split2)
print_message_with_time("########### Split2 model training completed, strat computing velocity ###########")
diff_mat_split2 = compute_sctour_velocity(tnode_split2, timestep=1/100)
print_message_with_time("########### Split2 velocity computed, start writing data ###########")
adata_split2.layers['velocity'] = diff_mat_split2
adata_split2.write(data_folder+'v2_erythroid/sct/adata_seed317_split2.h5ad')
tnode_split2.save_model(save_dir=data_folder+'v2_erythroid/sct/', save_prefix='tnode_seed317_split2')
print_message_with_time("########### Split2 data wrote ###########")


