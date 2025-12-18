# # scvelo virtual environment

import scvelo as scv
import numpy as np
import scanpy as sc
import math
import matplotlib.pyplot as plt
import scipy.io
import anndata as ad
import pandas as pd
import datetime

import sys
sys.path.append('/home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/veloUncertainty')
from countsplit import *
from v4_functions import print_message_with_time

data_folder = '/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/simulation/'

def run_countsplit_with_overdispersion(S,U,split_seed,overdisp_S,overdisp_U):
    print_message_with_time("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisp_S)
    u1, u2  = countsplit(U,overdisps=overdisp_U)
    return [[s1,u1],[s2,u2]]

def create_adata(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.obs = pd.DataFrame(index=adata_total.obs.index)
    for obs_col in adata_total.obs.columns:
        adata_split.obs[obs_col] = adata_total.obs[obs_col].copy()
    adata_split.var = pd.DataFrame(index=adata_total.var.index)
    for var_col in adata_total.var.columns:
        adata_split.var[var_col] = adata_total.var[var_col].copy()
    adata_split.layers['true_velocity'] = adata_total.layers['true_velocity'].copy()
    
    sc.pp.pca(adata_split)
    sc.pp.neighbors(adata_split, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata_split, n_pcs=None, n_neighbors=None)
    return adata_split

def countsplit_and_create_adata(S,U,total,split_seed,overdisp_S,overdisp_U):
    print_message_with_time("########### Running the function for overdispersion estimation and countsplitting")
    split1,split2 = run_countsplit_with_overdispersion(S=S,U=U,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata(split1[0],split1[1],total)
    adata2 = create_adata(split2[0],split2[1],total)
    return adata1,adata2


stochastic = []  # still collect in a list
for seed in range(50):
    print(f"Seed: {seed}")
    
    adata = ad.read_h5ad(data_folder + f"/adata_seed{seed}.h5ad")
    
    print_message_with_time("########### Reading data")
    overdisp_S = np.array(pd.read_csv(data_folder+'overdisp_S_seed' + str(seed) + '.csv')['x'])
    overdisp_U = np.array(pd.read_csv(data_folder+'overdisp_U_seed' + str(seed) + '.csv')['x'])
    
    S_mat = adata.layers['spliced'].copy()
    U_mat = adata.layers['unspliced'].copy()
    
    split_seed = 12122025
    print_message_with_time("########### Starting countsplit")
    adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)
    
    scv.tl.recover_dynamics(adata_split1)
    scv.tl.velocity(adata_split1, mode='stochastic')
    scv.tl.velocity_graph(adata_split1)
    
    scv.tl.recover_dynamics(adata_split2)
    scv.tl.velocity(adata_split2, mode='stochastic')
    scv.tl.velocity_graph(adata_split2)
    
    v1 = adata_split1.layers['velocity']
    v2 = adata_split2.layers['velocity']
    
    # Ensure inputs are numpy arrays (in case they are passed as other formats)
    v1 = np.array(v1)
    v2 = np.array(v2)
    
    # Safety check
    if v1.shape != v2.shape:
        raise ValueError(f"Shape mismatch: v1 is {v1.shape} but v2 is {v2.shape}")
    
    # --- Compute Cosine Similarity ---
    # Dot product per cell
    dot_product = np.sum(v1 * v2, axis=1)
    
    # Magnitudes (L2 norm)
    norm1 = np.linalg.norm(v1, axis=1)
    norm2 = np.linalg.norm(v2, axis=1)
    
    # Calculate similarity with epsilon safety
    epsilon = 1e-8
    cosine_sims = dot_product / (norm1 * norm2 + epsilon)
    
    mean_sim = np.mean(cosine_sims)
    stochastic.append(mean_sim)

stochastic = np.array(stochastic) 

import pandas as pd
pd.Series(stochastic, name="stochastic").to_csv("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/simulation/50-trials_stochastic.csv", index=False)