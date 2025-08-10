import sctour as sct
import scanpy as sc
import numpy as np
import torch
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *
from v4_functions import print_message_with_time

def run_sct_ery_Mark(gene_set_name,split_version,split_seed,sct_seed=615):
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method_prefix = 'sct'
    method = method_prefix + '_' + gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    print_message_with_time("################## Read data")
    raw = read_raw_adata(dataset_short)
    adata = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_'+split_version+'.h5ad')
    adata.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
    adata.obsm['X_pcaOriginal'] = raw.obsm['X_pca'].copy()
    print_message_with_time("########### start to train model for "+split_version+' ')
    tnode = sct_train_and_return_tnode(adata, sct_seed)
    print_message_with_time("########### start to compute velocity for "+split_version+' ')
    print_message_with_time("########### "+split_version+' velocity computed, start to write data')
    adata.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+split_version+'.h5ad')
    tnode.save_model(save_dir=savedata_folder, save_prefix='tnode_'+dataset_short+'_'+method+'_'+split_version)
    print_message_with_time("########### "+split_version+' data wrote')

for i in range(5):
    split_seed = [317, 320, 323, 326, 329][i]
    grid_seed = [227, 230, 233, 236, 239][i]
    gene_set_name = 'nMark' + str(grid_seed)
    run_sct_ery_Mark(gene_set_name=gene_set_name,split_version='split1',split_seed=split_seed,sct_seed=615)
    run_sct_ery_Mark(gene_set_name=gene_set_name,split_version='split2',split_seed=split_seed,sct_seed=615)
    run_sct_ery_Mark(gene_set_name=gene_set_name,split_version='total',split_seed=split_seed,sct_seed=615)
    