import scanpy as sc

dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'scv'
split_seed = 317
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

adata1 = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes_shrunk_disp.h5ad')
adata2 = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes_shrunk_disp.h5ad')
adata = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/'+dataset_short+'_total_allgenes_shrunk_disp.h5ad')

## run scv
import scvelo as scv
import bbknn
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'/'+method+'/'

"""
def scv_sto_compute_velocity_ery(adata):
    import bbknn
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    ### batch correction
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("Batch correction done!")
    del sys.modules['bbknn']
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40) 
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata,n_jobs=8)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata,n_jobs=8)
"""

print_message_with_time('################## Run model on total')
adata.layers['spliced_original'] = adata.layers['spliced'].copy()
adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
scv_compute_velocity_ery(adata) 

print_message_with_time('################## Read model on split1')
adata1.layers['spliced_original'] = adata1.layers['spliced'].copy()
adata1.layers['unspliced_original'] = adata1.layers['unspliced'].copy()
scv_compute_velocity_ery(adata1)

print_message_with_time('################## Read model on split2')
adata2.layers['spliced_original'] = adata2.layers['spliced'].copy()
adata2.layers['unspliced_original'] = adata2.layers['unspliced'].copy()
scv_compute_velocity_ery(adata2)

raw = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad')
adata.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
adata1.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
adata2.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()

# write data
print_message_with_time('################## Write data')
adata.write_h5ad(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_shrunk_disp.h5ad')
adata1.write_h5ad(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_shrunk_disp.h5ad')
adata2.write_h5ad(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_shrunk_disp.h5ad')



######
dataset_long = 'erythroid'
dataset_short = 'ery'
split_seed = 317
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'/'+method+'/'

import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

celltype_label = get_celltype_label(dataset_short)

scv.tl.velocity_confidence(adata)
scv.tl.velocity_confidence(adata1)
scv.tl.velocity_confidence(adata2)

## velocity
plot_velocity_scv_utv(adata_in=adata,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label='celltype')
plot_velocity_scv_utv(adata_in=adata1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label='celltype')
plot_velocity_scv_utv(adata_in=adata2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label='celltype')

## cosine similarity
plot_cosine_similarity(adata_split1=adata1,adata_split2=adata2,adata_total=adata,dataset=dataset_short,
						method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(adata1, adata2, adata, dataset_short, method, fig_folder, split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=adata1, adata_split2=adata2, adata_total=adata, dataset=dataset_short,
										method=method, fig_folder=fig_folder, split_seed=split_seed, celltype_label=celltype_label)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=adata1, adata_split2=adata2, adata_total=adata, dataset=dataset_short,
											method=method, fig_folder=fig_folder, split_seed=split_seed, celltype_label=celltype_label)

## confidence
plot_veloConf_and_cosSim(adata_total=adata, adata_split1=adata1, adata_split2=adata2, dataset=dataset_short,
							method=method, fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(adata, dataset_short, method, fig_folder, split_seed)
plot_velo_conf_boxplot_by_celltype(adata, dataset_short, method, fig_folder, split_seed, celltype_label=None)




