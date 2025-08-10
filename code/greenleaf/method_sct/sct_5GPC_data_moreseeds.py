method = 'sct_GPC'
dataset_long = 'greenleaf'
dataset_short = 'glf'
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
def read_data_and_run_sct_GPC(adata, dataset_short,dataset_long,method,savedata_folder,split_version,split_seed,sct_seed):
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


for split_seed in[320,323,326,329]:
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    split1 = sc.read(data_folder+'seed'+str(split_seed)+'_glf_split1_GPC.h5ad')
    split2 = sc.read(data_folder+'seed'+str(split_seed)+'_glf_split2_GPC.h5ad')
    split_version = 'split1'
    read_data_and_run_sct_GPC(split1,dataset_short,dataset_long,method,savedata_folder,split_version,split_seed,sct_seed)
    split_version = 'split2'
    read_data_and_run_sct_GPC(split2,dataset_short,dataset_long,method,savedata_folder,split_version,split_seed,sct_seed)
    print('################# seed'+str(split_seed)+' done!')