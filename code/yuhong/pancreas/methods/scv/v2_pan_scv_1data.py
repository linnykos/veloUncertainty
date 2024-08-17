import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime
from sklearn.metrics.pairwise import cosine_similarity

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

dataset_long = "pancreas"
dataset_short = "pan"
method = "scv"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

print_message_with_time("################## Read data")
total = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
adata_split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/seed317_split1_allgenes.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/seed317_split2_allgenes.h5ad')
gene_names = total.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()
S_mat_total = total.layers['spliced'].copy()
U_mat_total = total.layers['unspliced'].copy()

## run model
def scv_compute_velocity_pancreas(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # used to be n_neighbors=10, saved in 'Nbr10/'
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

print_message_with_time("################## Run model on total")
scv_compute_velocity_pancreas(total)
positions_total = [positions_dict[gene] for gene in total.var.index]
total.layers['spliced_original'] = S_mat_total[:,positions_total]
total.layers['unspliced_original'] = U_mat_total[:,positions_total]

print_message_with_time("################## Read model on split1")
scv_compute_velocity_pancreas(adata_split1)
positions_split1 = [positions_dict[gene] for gene in adata_split1.var.index]
adata_split1.layers['spliced_original'] = S_mat_split1[:,positions_split1] 
adata_split1.layers['unspliced_original'] = U_mat_split1[:,positions_split1]

print_message_with_time("################## Read model on split2")
scv_compute_velocity_pancreas(adata_split2)
positions_split2 = [positions_dict[gene] for gene in adata_split2.var.index]
adata_split2.layers['spliced_original'] = S_mat_split2[:,positions_split2]
adata_split2.layers['unspliced_original'] = U_mat_split2[:,positions_split2]

# write data
print_message_with_time("################## Write data")
total.write_h5ad(data_folder+'v2_'+dataset_long+'/scv/adata_'+dataset_short+'_scv_total_v2.h5ad')
adata_split1.write_h5ad(data_folder+'v2_'+dataset_long+'/scv/adata_'+dataset_short+'_scv_split1_v2.h5ad')
adata_split2.write_h5ad(data_folder+'v2_'+dataset_long+'/scv/adata_'+dataset_short+'_scv_split2_v2.h5ad')

print_message_with_time("################## All done with the data")

import numpy as np
common_genes_filter = np.intersect1d(np.array(adata_split1.var.index), np.array(adata_split2.var.index))
print('Number of overlapped genes between splits = '+str(common_genes_filter.shape[0])) 
print('Number of overlapped genes between split1 and total = '+str(np.intersect1d(np.array(adata_split1.var.index), np.array(total.var.index)).shape[0])) # 589 -> 1335
print('Number of overlapped genes between split2 and total = '+str(np.intersect1d(np.array(adata_split2.var.index), np.array(total.var.index)).shape[0])) # 580 -> 1328

print("Number of genes used in velocity computation in split1 = "+str(np.sum(~np.isnan(adata_split1.layers['velocity'][0])))) 
print("Number of genes used in velocity computation in split2 = "+str(np.sum(~np.isnan(adata_split2.layers['velocity'][0])))) 
print("Number of genes used in velocity computation in total = "+str(np.sum(~np.isnan(total.layers['velocity'][0])))) 

velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0]))



