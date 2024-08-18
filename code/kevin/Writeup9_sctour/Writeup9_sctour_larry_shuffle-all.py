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

tnode = sct.train.Trainer(adata, 
                          loss_mode='nb', 
                          alpha_recon_lec=0.5, 
                          alpha_recon_lode=0.5)
tnode.train()
adata.obs['ptime'] = tnode.get_time()
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, 
                                         alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, 
                                            adata.obsm['X_TNODE'])

adata.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup0_sctour_larry_shuffle-all.h5ad")
tnode.save_model(save_dir="/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/", 
                 save_prefix="Writeup9_sctour_larry_tnode_shuffle-all")
print_message_with_time("########### Total data wrote")

