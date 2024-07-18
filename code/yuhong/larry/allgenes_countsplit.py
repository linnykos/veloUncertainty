import scvelo as scv
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

split_seed = 317

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

adata = ad.read_h5ad(data_folder+"v2_larry/larry.h5ad") # n_obs × n_vars = 49302 × 23420

S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

print_message_with_time("########### Starting countsplit")

def run_countsplit_with_overdispersion(S,U,split_seed):
    print_message_with_time("########### Estimating overdispersion parameters")
    overdisps_S = estimate_overdisps(S)
    overdisps_U = estimate_overdisps(U)
    print_message_with_time("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisps_S)
    u1, u2  = countsplit(U,overdisps=overdisps_U)
    return [[s1,u1],[s2,u2]]

def create_adata_larry(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
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
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata_larry(split1[0],split1[1],total)
    adata2 = create_adata_larry(split2[0],split2[1],total)
    return adata1,adata2

adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=317)

adata_split1.write(data_folder+'v2_larry/larry_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v2_larry/larry_split2_allgenes.h5ad')

total = adata.copy()
del total.obsm
del total.layers['ambiguous']
del total.layers['matrix']
total.write(data_folder+'v2_larry/larry_total_allgenes.h5ad')