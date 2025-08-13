import scanpy as sc
import pandas as pd
import scvelo as scv
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

method_prefix = 'scv'
gene_set_name = 'Mark'
method = method_prefix + '_' + gene_set_name

for split_seed in [317,320,323,326,329]:
    total = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_total_allgenes.h5ad')
    split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_split1_allgenes.h5ad')
    split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_split2_allgenes.h5ad')
    scv_compute_velocity_ery(total) 
    scv_compute_velocity_ery(split1) 
    scv_compute_velocity_ery(split2) 
    total.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')

print('####################### All done.')
