dataset_long = 'erythroid'
dataset_short = 'ery'
method_prefix = 'scv'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

import scanpy as sc
import pandas as pd
import scvelo as scv
import numpy as np
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

##### process data
for i in range(4):
    split_seed = [320, 323, 326, 329][i]
    grid_seed = 227
    gene_set_name = 'nMark' + str(grid_seed)
    method = method_prefix + '_' + gene_set_name
    total = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_total.h5ad')
    ## manually adjust total object
    total.obsm['X_umapOriginal'] = total.obsm['X_umap'].copy()
    total.obsm['X_pcaOriginal'] = total.obsm['X_pca'].copy()
    del total.obsm['X_umap']
    del total.obsm['X_pca']
    split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_split1.h5ad')
    split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_split2.h5ad')
    scv_compute_velocity_ery(total) 
    scv_compute_velocity_ery(split1) 
    scv_compute_velocity_ery(split2) 
    total.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
    print('##################### split_seed'+str(split_seed)+' data done!')


