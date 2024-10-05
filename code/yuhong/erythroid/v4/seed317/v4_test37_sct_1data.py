split_seed = 317
method = 'sct'
dataset_long = 'erythroid'
dataset_short = 'ery'
sct_seed = 615

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_test/v4_'+dataset_long+'/'
savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'

import sctour as sct
import scanpy as sc
import numpy as np
import torch
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *

######################################################
import random
def sct_train_and_return_tnode(adata, sct_seed=615):
    import sctour as sct
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    adata.X = adata.X.astype(np.float32)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)
    tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
    return tnode


raw = read_raw_adata(dataset_short)

split_version = '37-3'
adata = sc.read_h5ad(data_folder+'seed317_ery_allgenes_'+split_version+'.h5ad')
if not 'larry' in dataset_long:
    adata.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
    adata.obsm['X_pcaOriginal'] = raw.obsm['X_pca'].copy()
print_message_with_time("########### start to train model for "+split_version+' ')
tnode = sct_train_and_return_tnode(adata, sct_seed)
print_message_with_time("########### start to compute velocity for "+split_version+' ')
adata.layers['spliced_original'] = adata.layers['spliced'].copy()
adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
adata.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+split_version+'_v4.h5ad')
tnode.save_model(save_dir=savedata_folder, save_prefix='tnode_'+dataset_short+'_'+method+'_'+split_version+'_v4')
print_message_with_time("########### "+split_version+' data wrote')

split_version = '37-7'
adata = sc.read_h5ad(data_folder+'seed317_ery_allgenes_'+split_version+'.h5ad')
if not 'larry' in dataset_long:
    adata.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
    adata.obsm['X_pcaOriginal'] = raw.obsm['X_pca'].copy()
print_message_with_time("########### start to train model for "+split_version+' ')
tnode = sct_train_and_return_tnode(adata, sct_seed)
print_message_with_time("########### start to compute velocity for "+split_version+' ')
adata.layers['spliced_original'] = adata.layers['spliced'].copy()
adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
adata.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+split_version+'_v4.h5ad')
tnode.save_model(save_dir=savedata_folder, save_prefix='tnode_'+dataset_short+'_'+method+'_'+split_version+'_v4')
print_message_with_time("########### "+split_version+' data wrote')

