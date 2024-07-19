import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from sctour_misc import *

method = 'sct'
dataset_long = 'larry'
dataset_short = 'larry'

sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

######################################################

print_message_with_time("################## Read data")
#raw = sc.read_h5ad(data_folder+"v2_larry/larry.h5ad")
adata_split2 = sc.read_h5ad(data_folder+'v2_larry/larry_split2_allgenes.h5ad')
gene_names = adata_split2.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()

######################################################

def sct_train_and_return_tnode(adata, sct_seed=615):
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

print_message_with_time("########### Start to train model for split2")
tnode_split2 = sct_train_and_return_tnode(adata_split2)
print_message_with_time("########### Start to compute velocity for split2")
diff_mat_split2 = compute_sctour_velocity(tnode_split2, timestep=1/100)
print_message_with_time("########### Split2 velocity computed, start to write data")
adata_split2.layers['velocity'] = diff_mat_split2
adata_split2.layers['spliced_original'] = adata_split2.layers['spliced'].copy()
adata_split2.layers['unspliced_original'] = adata_split2.layers['unspliced'].copy()
adata_split2.write(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')
tnode_split2.save_model(save_dir=data_folder+'v2_'+dataset_long+'/'+method+'/', save_prefix='tnode_'+dataset_short+'_'+method+'_split2_v2')
print_message_with_time("########### Split2 data wrote")


