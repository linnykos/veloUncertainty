# unitvelo virtual environment

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import datetime
import random
import numpy as np
from scipy.sparse import csr_matrix

sct_seed = 615

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1

os.environ["TF_USE_LEGACY_KERAS"]="1"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

adata = sc.read_h5ad("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry.h5ad")

random.seed(sct_seed)
np.random.seed(sct_seed)

# Check if the layers exist before trying to delete them
layers_to_remove = ['ambiguous', 'matrix']

for layer in layers_to_remove:
    if layer in adata.layers:
        del adata.layers[layer]

# Verify the layers have been removed
print(adata.layers)

################

sc.pp.highly_variable_genes(adata, 
                            flavor='seurat_v3', 
                            n_top_genes=2000, 
                            subset=True)

################

# shuffle the genes 

# Ensure that the data matrices are dense if they are sparse
adata_X_dense = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
adata_spliced_dense = adata.layers['spliced'].toarray() if not isinstance(adata.layers['spliced'], np.ndarray) else adata.layers['spliced']
adata_unspliced_dense = adata.layers['unspliced'].toarray() if not isinstance(adata.layers['unspliced'], np.ndarray) else adata.layers['unspliced']

# Shuffle values independently for each gene (column)
for i in range(adata_X_dense.shape[1]):
    shuffle_order = np.random.permutation(adata_X_dense.shape[0])  # Generate a random shuffle order for the rows (cells)
    
    # Apply the shuffle order to adata.X, adata.layers['spliced'], and adata.layers['unspliced']
    adata_X_dense[:, i] = adata_X_dense[shuffle_order, i]
    adata_spliced_dense[:, i] = adata_spliced_dense[shuffle_order, i]
    adata_unspliced_dense[:, i] = adata_unspliced_dense[shuffle_order, i]

# Replace the original data with the shuffled data
adata.X = csr_matrix(adata_X_dense)
adata.layers['spliced'] = csr_matrix(adata_spliced_dense)
adata.layers['unspliced'] = csr_matrix(adata_unspliced_dense)

################
adata = utv.run_model(adata,
                      'state_info', 
                      config_file=velo_config)

adata.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup9_unitvelo_larry_shuffle-all.h5ad")
print_message_with_time("########### Total data wrote")

