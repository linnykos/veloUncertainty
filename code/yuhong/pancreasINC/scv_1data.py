import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
import datetime
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

dataset_long = "pancreasINC"
dataset_short = "panINC"
method = "scv"

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_pancreasINC/'+method+'/'

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

# pancreasINC_split1_allgenes.h5ad  pancreasINC_split2_allgenes.h5ad  pancreasINC_total_allgenes.h5ad

print_message_with_time("################## Read data")
#raw = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
total = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_total_allgenes.h5ad')
adata_split1 = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_split1_allgenes.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_split2_allgenes.h5ad')
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
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # used to be n_neighbors=10
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

print_message_with_time("################## Run model on total")
scv_compute_velocity_pancreas(total) # 18:12-18:22
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

### common genes being filtered out
common_genes = np.intersect1d(np.array(adata_split1.var.index), np.array(adata_split2.var.index)) 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 1442
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(adata_split1.var.index),np.array(total.var.index)).shape[0])) # 1248
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(adata_split2.var.index),np.array(total.var.index)).shape[0])) # 1269

plot_gene_correlation_between_splits(adata1=adata_split1,adata2=adata_split2,fig_folder=fig_folder,fig_path="split_cor/"+dataset_short+"_"+method+"_corr_between_splits_colorSame.png")
