import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

random.seed(sct_seed)
np.random.seed(sct_seed)

################

# Check if the layers exist before trying to delete them
layers_to_remove = ['ambiguous', 'matrix']

for layer in layers_to_remove:
    if layer in adata.layers:
        del adata.layers[layer]

# Verify the layers have been removed
print(adata.layers)

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

S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

overdisps_S = estimate_overdisps(S_mat)
overdisps_U = estimate_overdisps(U_mat)

np.random.seed(sct_seed)
S_split1, S_split2  = countsplit(S_mat,overdisps=overdisps_S)
np.random.seed(sct_seed)
U_split1, U_split2  = countsplit(U_mat,overdisps=overdisps_U)

adata_split1 = ad.AnnData(X=S_split1.astype(np.float32))
adata_split1.layers["spliced"] = S_split1
adata_split1.layers["unspliced"] = U_split1
adata_split1.obs = pd.DataFrame(index=adata.obs.index)
for obs_col in adata.obs.columns:
    adata_split1.obs[obs_col] = adata.obs[obs_col].copy()

adata_split1.var = pd.DataFrame(index=adata.var.index)
for var_col in adata.var.columns:
    adata_split1.var[var_col] = adata.var[var_col].copy()

###

adata_split2 = ad.AnnData(X=S_split2.astype(np.float32))
adata_split2.layers["spliced"] = S_split2
adata_split2.layers["unspliced"] = U_split2
adata_split2.obs = pd.DataFrame(index=adata.obs.index)
for obs_col in adata.obs.columns:
    adata_split2.obs[obs_col] = adata.obs[obs_col].copy()

adata_split2.var = pd.DataFrame(index=adata.var.index)
for var_col in adata.var.columns:
    adata_split2.var[var_col] = adata.var[var_col].copy()

adata_split1.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup9_larry_split1.h5ad")
adata_split2.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup9_larry_split2.h5ad")
print_message_with_time("########### Finish splitting")