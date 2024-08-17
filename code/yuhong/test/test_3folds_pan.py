import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from v2_functions import print_message_with_time

dataset_short = 'pan'
dataset_long = 'pancreas'
split_seed = 317


print_message_with_time("########### Starting countsplit")

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

adata = sc.read(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

def run_countsplit_with_overdispersion_3folds(S,U,split_seed):
    print_message_with_time("########### Estimating overdispersion parameters")
    overdisps_S = estimate_overdisps(S)
    overdisps_U = estimate_overdisps(U)
    print_message_with_time("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2, s3  = countsplit(S,overdisps=overdisps_S,folds=3)
    u1, u2, u3 = countsplit(U,overdisps=overdisps_U,folds=3)
    return [[s1,u1],[s2,u2],[s3,u3]]

def create_adata_pancreas(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.uns = {'clusters_colors':adata.uns['clusters_colors'].copy()}
    adata_split.obs = pd.DataFrame(index=adata_total.obs.index)
    for obs_col in adata_total.obs.columns:
        adata_split.obs[obs_col] = adata_total.obs[obs_col].copy()
    adata_split.var = pd.DataFrame(index=adata_total.var.index)
    for var_col in adata_total.var.columns:
        adata_split.var[var_col] = adata_total.var[var_col].copy()
    return adata_split

def countsplit_and_create_adata(S,U,total,split_seed):
    print_message_with_time("########### Running the function for overdispersion estimation and countsplitting")
    split1,split2,split3 = run_countsplit_with_overdispersion_3folds(S=S,U=U,split_seed=split_seed)
    # split1
    counts_adata1 = ad.AnnData(X=split1[0].astype(np.float32))
    counts_adata1.layers["spliced"] = split1[0]
    counts_adata1.layers["unspliced"] = split1[1]
    # split2
    counts_adata2 = ad.AnnData(X=split2[0].astype(np.float32))
    counts_adata2.layers["spliced"] = split2[0]
    counts_adata2.layers["unspliced"] = split2[1]
    # split3
    counts_adata3 = ad.AnnData(X=split3[0].astype(np.float32))
    counts_adata3.layers["spliced"] = split3[0]
    counts_adata3.layers["unspliced"] = split3[1]
    print_message_with_time("########### Creating split adata objects")
    adata1 = create_adata_pancreas(split1[0],split1[1],total)
    adata2 = create_adata_pancreas(split2[0],split2[1],total)
    adata3 = create_adata_pancreas(split3[0],split3[1],total)
    return adata1,adata2,adata3

adata_split1,adata_split2,adata_split3 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=adata,split_seed=317)

adata_split1.write(data_folder+'test/'+dataset_short+'_split1_3folds_allgenes.h5ad')
adata_split2.write(data_folder+'test/'+dataset_short+'_split2_3folds_allgenes.h5ad')
adata_split3.write(data_folder+'test/'+dataset_short+'_split3_3folds_allgenes.h5ad')
