split_seed = 320

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

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

def run_countsplit_with_overdispersion(S,U,split_seed,overdisp_S,overdisp_U):
    print_message_with_time("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisp_S)
    u1, u2  = countsplit(U,overdisps=overdisp_U)
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

def countsplit_and_create_adata(S,U,total,split_seed,overdisp_S,overdisp_U):
    print_message_with_time("########### Running the function for overdispersion estimation and countsplitting")
    split1,split2 = run_countsplit_with_overdispersion(S=S,U=U,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata_larry(split1[0],split1[1],total)
    adata2 = create_adata_larry(split2[0],split2[1],total)
    return adata1,adata2

print_message_with_time("########### Reading data")
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

overdisp_S = np.array(pd.read_csv(data_folder+'v4_larry/larry_overdisp_S.csv')['x'])
overdisp_U = np.array(pd.read_csv(data_folder+'v4_larry/larry_overdisp_U.csv')['x'])

adata = ad.read_h5ad(data_folder+"v4_larry/larry.h5ad") # n_obs × n_vars = 49302 × 23420
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("########### Starting countsplit")
adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)

adata_split1.write(data_folder+'v4_larry/seed'+str(split_seed)+'_larry_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v4_larry/seed'+str(split_seed)+'_larry_split2_allgenes.h5ad')

print_message_with_time("########### All done")
