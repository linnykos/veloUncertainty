import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

dataset_long = "pancreas"
dataset_short = "pan"
method = "scv"

print_message_with_time("################## Read data")
total = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
adata_split1 = sc.read_h5ad(data_folder+'test/'+dataset_short+'_split1_3folds_allgenes.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'test/'+dataset_short+'_split2_3folds_allgenes.h5ad')
adata_split3 = sc.read_h5ad(data_folder+'test/'+dataset_short+'_split3_3folds_allgenes.h5ad')

gene_names = total.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()
S_mat_split3 = adata_split3.layers['spliced'].copy()
U_mat_split3 = adata_split3.layers['unspliced'].copy()

S_mat_total = total.layers['spliced'].copy()
U_mat_total = total.layers['unspliced'].copy()

## run model
def scv_compute_velocity_pancreas(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) 
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

print_message_with_time("################## Run model on total")
scv_compute_velocity_pancreas(total)
positions_total = [positions_dict[gene] for gene in total.var.index]
total.layers['spliced_original'] = S_mat_total[:,positions_total]
total.layers['unspliced_original'] = U_mat_total[:,positions_total]

print_message_with_time("################## Run model on split1")
scv_compute_velocity_pancreas(adata_split1)
positions_split1 = [positions_dict[gene] for gene in adata_split1.var.index]
adata_split1.layers['spliced_original'] = S_mat_split1[:,positions_split1] 
adata_split1.layers['unspliced_original'] = U_mat_split1[:,positions_split1]

print_message_with_time("################## Run model on split2")
scv_compute_velocity_pancreas(adata_split2)
positions_split2 = [positions_dict[gene] for gene in adata_split2.var.index]
adata_split2.layers['spliced_original'] = S_mat_split2[:,positions_split2]
adata_split2.layers['unspliced_original'] = U_mat_split2[:,positions_split2]

print_message_with_time("################## Run model on split3")
scv_compute_velocity_pancreas(adata_split3)
positions_split3 = [positions_dict[gene] for gene in adata_split3.var.index]
adata_split3.layers['spliced_original'] = S_mat_split3[:,positions_split3]
adata_split3.layers['unspliced_original'] = U_mat_split3[:,positions_split3]


# write data
print_message_with_time("################## Write data")
total.write_h5ad(data_folder+'test/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
adata_split1.write_h5ad(data_folder+'test/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
adata_split2.write_h5ad(data_folder+'test/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
adata_split3.write_h5ad(data_folder+'test/'+method+'/adata_'+dataset_short+'_'+method+'_split3.h5ad')

print_message_with_time("################## All done with the data")



