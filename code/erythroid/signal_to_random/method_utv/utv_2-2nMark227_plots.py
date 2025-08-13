import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import plot_velocity_scv_utv
from v4_functions import *



def compute_umap_utv_ery_nMark(gene_set_name, split_seed):
    dataset_long = 'erythroid'
    dataset_short = 'ery'
    method_prefix = 'utv'
    method = method_prefix + '_' + gene_set_name
    celltype_label = get_celltype_label(dataset_short)
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    total.obsm['X_umapOriginal'] = total.obsm['X_umap'].copy()
    total.obsm['X_pcaOriginal'] = total.obsm['X_pca'].copy()
    del total.obsm['X_umap']
    del total.obsm['X_pca']
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
    ## compute umap
    compute_umap_ery(total)
    compute_umap_ery(split1)
    compute_umap_ery(split2)
    scv.tl.velocity_graph(total)
    scv.tl.velocity_graph(split1)
    scv.tl.velocity_graph(split2)
    total.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_outputAdded.h5ad')
    split1,write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded.h5ad')
    split2,write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded.h5ad')

split_seed = 317
grid_seed = 227
gene_set_name = 'nMark' + str(grid_seed)
compute_umap_utv_ery_nMark(gene_set_name=gene_set_name, split_seed=split_seed)
print('##################### seed317 done')

split_seed = 320
grid_seed = 230
gene_set_name = 'nMark' + str(grid_seed)
compute_umap_utv_ery_nMark(gene_set_name='nMark', split_seed=320)
print('##################### seed320 done')

split_seed = 323
grid_seed = 233
gene_set_name = 'nMark' + str(grid_seed)
compute_umap_utv_ery_nMark(gene_set_name='nMark', split_seed=323)
print('##################### seed323 done')

split_seed = 326
grid_seed = 236
gene_set_name = 'nMark' + str(grid_seed)
compute_umap_utv_ery_nMark(gene_set_name='nMark', split_seed=326)
print('##################### seed326 done')

split_seed = 329
grid_seed = 239
gene_set_name = 'nMark' + str(grid_seed)
compute_umap_utv_ery_nMark(gene_set_name='nMark', split_seed=329)
print('##################### seed329 done')

