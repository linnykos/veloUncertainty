import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

split_seed = 317

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

print_message_with_time("########### Starting countsplit")

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

adata = sc.read(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

def run_countsplit_with_overdispersion(S,U,split_seed):
    print_message_with_time("########### Estimating overdispersion parameters")
    overdisps_S = estimate_overdisps(S)
    overdisps_U = estimate_overdisps(U)
    print_message_with_time("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisps_S)
    u1, u2  = countsplit(U,overdisps=overdisps_U)
    return [[s1,u1],[s2,u2]]#[s1,s2,u1,u2]

def create_adata_pancreas(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.obs = pd.DataFrame(index=adata_total.obs.index)
    adata_split.obs['clusters'] = adata_total.obs['clusters'].copy()
    adata_split.var = pd.DataFrame(index=adata_total.var.index)
    adata_split.var['high_variable_genes'] = adata_total.var['high_variable_genes'].index.copy()
    adata_split.uns = {'clusters_colors':adata.uns['clusters_colors'].copy()}
    adata_split.obsm['X_pcaOriginal'] = adata_total.obsm['X_pca'].copy()
    adata_split.obsm['X_umapOriginal'] = adata_total.obsm['X_umap'].copy()
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
    counts_adata1.write(data_folder+'v2_pancreas/counts_seed317_split1_allgenes.h5ad')
    counts_adata2.write(data_folder+'v2_pancreas/counts_seed317_split2_allgenes.h5ad')
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata_pancreas(split1[0],split1[1],total)
    adata2 = create_adata_pancreas(split2[0],split2[1],total)
    return adata1,adata2

adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=317)

adata_split1.write(data_folder+'v2_pancreas/seed317_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v2_pancreas/seed317_split2_allgenes.h5ad')
