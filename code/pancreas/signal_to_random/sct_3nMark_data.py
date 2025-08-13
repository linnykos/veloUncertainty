dataset_long = 'pancreas'
dataset_short = 'pan'
sct_seed = 615
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

import sctour as sct
import scanpy as sc
import numpy as np
import torch
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *

######################################################
def read_data_and_run_sct_nMark(adata, dataset_short,dataset_long,method,savedata_folder,split_version,sct_seed):
    print_message_with_time("################## Read data")
    raw = read_raw_adata(dataset_short)
    if (not 'larry' in dataset_long) and (not dataset_long=='greenleaf'):
        adata.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
        adata.obsm['X_pcaOriginal'] = raw.obsm['X_pca'].copy()
    print_message_with_time("########### start to train model for "+split_version+' ')
    tnode = sct_train_and_return_tnode(adata, sct_seed)
    print_message_with_time("########### start to compute velocity for "+split_version+' ')
    print_message_with_time("########### "+split_version+' velocity computed, start to write data')
    adata.layers['spliced_original'] = adata.layers['spliced'].copy()
    adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
    adata.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+split_version+'_v4.h5ad')
    tnode.save_model(save_dir=savedata_folder, save_prefix='tnode_'+dataset_short+'_'+method+'_'+split_version+'_v4')
    print_message_with_time("########### "+split_version+' data wrote')


for i in range(5):
    split_seed = [317, 320, 323, 326, 329][i]
    grid_seed = [227, 230, 233, 236, 239][i]
    gene_set_name = 'nMark'+str(grid_seed)
    method = 'sct_'+gene_set_name
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    total = sc.read(data_folder+dataset_short+'_total_'+gene_set_name+'.h5ad')
    split1 = sc.read(data_folder+dataset_short+'_split1_'+gene_set_name+'.h5ad')
    split2 = sc.read(data_folder+dataset_short+'_split2_'+gene_set_name+'.h5ad')
    read_data_and_run_sct_nMark(split1,dataset_short,dataset_long,method,savedata_folder,split_version='split1',sct_seed=sct_seed)
    read_data_and_run_sct_nMark(split2,dataset_short,dataset_long,method,savedata_folder,split_version='split2',sct_seed=sct_seed)
    read_data_and_run_sct_nMark(total,dataset_short,dataset_long,method,savedata_folder,split_version='total',sct_seed=sct_seed)
    print('################### seed'+str(split_seed)+' done')


