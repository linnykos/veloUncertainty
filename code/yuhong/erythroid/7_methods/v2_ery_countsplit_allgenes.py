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

adata = sc.read(data_folder+"Gastrulation/erythroid_lineage.h5ad")
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

def create_adata_erythroid(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.obs = pd.DataFrame(index=adata_total.obs.index)
    adata_split.obs['celltype'] = adata_total.obs['celltype'].copy()
    adata_split.obs['sequencing.batch'] = adata_total.obs['sequencing.batch'].copy()
    adata_split.var = pd.DataFrame(index=adata_total.var.index)
    adata_split.var['Accession'] = adata_total.var['Accession'].index.copy()
    adata_split.uns['celltype_colors'] = {'celltype_colors':adata.uns['celltype_colors'].copy()}
    adata_split.obsm['X_pcaOriginal'] = adata_total.obsm['X_pca'].copy()
    adata_split.obsm['X_umapOriginal'] = adata_total.obsm['X_umap'].copy()
    return adata_split

def countsplit_and_create_adata(S,U,total,split_seed):
    print_message_with_time("########### Running the function for overdispersion estimation and countsplitting")
    split1,split2 = run_countsplit_with_overdispersion(S=S,U=U,split_seed=split_seed)
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata_erythroid(split1[0],split1[1],total)
    adata2 = create_adata_erythroid(split2[0],split2[1],total)
    return adata1,adata2

adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=317)

adata_split1.write(data_folder+'v2_erythroid/seed317_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v2_erythroid/seed317_split2_allgenes.h5ad')
