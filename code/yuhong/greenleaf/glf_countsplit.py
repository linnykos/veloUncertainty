dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

import anndata as ad
adata = ad.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 

import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from v4_functions import print_message_with_time

def run_countsplit_with_overdispersion(S,U,split_seed,overdisp_S,overdisp_U):
    print_message_with_time("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisp_S)
    u1, u2  = countsplit(U,overdisps=overdisp_U)
    return [[s1,u1],[s2,u2]]

def create_adata_greenleaf(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.uns = {'cluster_name_colors':adata.uns['cluster_name_colors'].copy()}
    #adata_split.obsm['X_pcaOriginal'] = adata_total.obsm['X_pca'].copy()
    adata_split.obsm['X_umapOriginal'] = adata_total.obsm['X_umap_greenleaf'].copy()
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
    adata1 = create_adata_greenleaf(split1[0],split1[1],total)
    adata2 = create_adata_greenleaf(split2[0],split2[1],total)
    return adata1,adata2

print_message_with_time("########### Reading data")
overdisp_S = np.array(pd.read_csv(data_folder+'v4_'+dataset_long+'/'+'greenleaf_overdisp_S.csv')['x'])
overdisp_U = np.array(pd.read_csv(data_folder+'v4_'+dataset_long+'/'+'greenleaf_overdisp_U.csv')['x'])

S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

## 317
split_seed = 317
print_message_with_time("########### Starting countsplit "+str(split_seed))
adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)

adata_split1.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

## 320
split_seed = 320
print_message_with_time("########### Starting countsplit "+str(split_seed))
adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)

adata_split1.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

## 323
split_seed = 323
print_message_with_time("########### Starting countsplit "+str(split_seed))
adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)

adata_split1.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

## 326
split_seed = 326
print_message_with_time("########### Starting countsplit "+str(split_seed))
adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)

adata_split1.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

## 329
split_seed = 329
print_message_with_time("########### Starting countsplit "+str(split_seed))
adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)

adata_split1.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

print_message_with_time("########### All done")
