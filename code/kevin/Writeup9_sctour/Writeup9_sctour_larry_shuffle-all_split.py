import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

adata = sc.read_h5ad("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry.h5ad")

torch.manual_seed(sct_seed)
random.seed(sct_seed)
np.random.seed(sct_seed)

################

# Check if the layers exist before trying to delete them
layers_to_remove = ['spliced', 'unspliced', 'ambiguous', 'matrix']

for layer in layers_to_remove:
    if layer in adata.layers:
        del adata.layers[layer]

# Verify the layers have been removed
print(adata.layers)

################

# shuffle the genes 

# Assuming adata.X is in a sparse matrix format, convert it to a dense matrix
adata.X = adata.X.toarray()  # Convert to dense matrix if it's sparse

# Shuffle values independently for each gene (column)
for i in range(adata.X.shape[1]):
    np.random.shuffle(adata.X[:, i])

# If you want to convert back to a sparse matrix, you can do so (optional)
adata.X = csr_matrix(adata.X)

################

adata.X = adata.X.astype(np.float32)
sc.pp.calculate_qc_metrics(adata, 
                           percent_top=None, 
                           log1p=False, 
                           inplace=True)
sc.pp.highly_variable_genes(adata, 
                            flavor='seurat_v3', 
                            n_top_genes=2000, 
                            subset=True)

################

# now split
# code from: https://github.com/linnykos/veloUncertainty/blob/main/code/yuhong/larry/allgenes_countsplit.py

X_mat = adata.X.copy()
overdisps = estimate_overdisps(X_mat)

np.random.seed(sct_seed)
X_split1, X_split2  = countsplit(X_mat)

adata_split1 = ad.AnnData(X=X_split1.astype(np.float32))
adata_split1.obs = pd.DataFrame(index=adata.obs.index)
for obs_col in adata.obs.columns:
    adata_split1.obs[obs_col] = adata.obs[obs_col].copy()

adata_split1.var = pd.DataFrame(index=adata.var.index)
for var_col in adata.var.columns:
    adata_split1.var[var_col] = adata.var[var_col].copy()

###

adata_split2 = ad.AnnData(X=X_split2.astype(np.float32))
adata_split2.obs = pd.DataFrame(index=adata.obs.index)
for obs_col in adata.obs.columns:
    adata_split2.obs[obs_col] = adata.obs[obs_col].copy()

adata_split2.var = pd.DataFrame(index=adata.var.index)
for var_col in adata.var.columns:
    adata_split2.var[var_col] = adata.var[var_col].copy()

################

# first train on split 1

tnode_split1 = sct.train.Trainer(adata_split1, 
                                 loss_mode='nb', 
                                 alpha_recon_lec=0.5, 
                                 alpha_recon_lode=0.5)
tnode_split1.train()
adata_split1.obs['ptime'] = tnode_split1.get_time()
mix_zs, zs, pred_zs = tnode_split1.get_latentsp(alpha_z=0.5, 
                                                alpha_predz=0.5)
adata_split1.obsm['X_TNODE'] = mix_zs
adata_split1.obsm['X_VF'] = tnode_split1.get_vector_field(adata_split1.obs['ptime'].values, 
                                                          adata_split1.obsm['X_TNODE'])

adata_split1.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup0_sctour_larry_shuffle-all_split1.h5ad")
tnode_split1.save_model(save_dir="/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/", 
                        save_prefix="Writeup9_sctour_larry_tnode_shuffle-all_split1")
print_message_with_time("########### Split 1 wrote")

# train on split 2

tnode_split2 = sct.train.Trainer(adata_split2, 
                                 loss_mode='nb', 
                                 alpha_recon_lec=0.5, 
                                 alpha_recon_lode=0.5)
tnode_split2.train()
adata_split2.obs['ptime'] = tnode_split2.get_time()
mix_zs, zs, pred_zs = tnode_split2.get_latentsp(alpha_z=0.5, 
                                                alpha_predz=0.5)
adata_split2.obsm['X_TNODE'] = mix_zs
adata_split2.obsm['X_VF'] = tnode_split2.get_vector_field(adata_split2.obs['ptime'].values, 
                                                          adata_split2.obsm['X_TNODE'])

adata_split2.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup0_sctour_larry_shuffle-all_split2.h5ad")
tnode_split2.save_model(save_dir="/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/", 
                        save_prefix="Writeup9_sctour_larry_tnode_shuffle-all_split2")
print_message_with_time("########### Split 2 wrote")