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

def addvelo_utv_ery_Mark(gene_set_name, split_seed, total_velo=False):
    dataset_long = 'erythroid'
    dataset_short = 'ery'
    method_prefix = 'utv'
    method = method_prefix + '_' + gene_set_name
    celltype_label = get_celltype_label(dataset_short)
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    if total_velo: total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
    ## compute umap
    if total_velo: compute_umap_ery(total)
    compute_umap_ery(split1)
    compute_umap_ery(split2)
    if total_velo: scv.tl.velocity_graph(total)
    scv.tl.velocity_graph(split1)
    scv.tl.velocity_graph(split2)
    ## write data
    if total_velo: total.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_outputAdded.h5ad')
    split1.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded.h5ad')
    split2.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded.h5ad')
    print('###################### seed'+str(split_seed)+' done.')
    
total_velo = True
for i in range(5):
    split_seed = [317,320,323,326,329][i]
    #split_seed = [323,326,329][i]
    addvelo_utv_ery_Mark(gene_set_name='Mark', split_seed=split_seed, total_velo=total_velo)
    total_velo = False

print('###################### All done.')
