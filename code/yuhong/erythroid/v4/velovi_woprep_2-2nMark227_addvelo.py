method_prefix = 'velovi_woprep'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_velovi import *
from v4_functions import *

def save_adata_outputAdded(data_version, method_prefix, gene_set_name, split_seed):
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method = method_prefix + '_' + gene_set_name
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'
    #fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    ## read data
    adata = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'.h5ad')
    vae = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'.pt', adata)
    if data_version=='total': 
        adata.obsm['X_umapOriginal'] = adata.obsm['X_umap']
        del adata.obsm['X_umap']
    ## add velovi outputs to adata
    print_message_with_time("############## Add velovi outputs to adata")
    add_velovi_outputs_to_adata(adata, vae)
    ######################################################
    ## compute umap
    print_message_with_time("############## Compute umap")
    compute_umap_ery(adata)
    ######################################################
    ## velocity graph
    print_message_with_time("############## Compute velocity graph")
    scv.tl.velocity_graph(adata)
    ######################################################
    ## save adata with velocity
    adata.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_outputAdded.h5ad')
    print('######################## output added!')

for i in range(5):
    split_seed = [317, 320, 323, 326, 329][i]
    grid_seed = [227, 230, 233, 236, 239][i]
    gene_set_name = 'nMark' + str(grid_seed)
    save_adata_outputAdded(data_version='total', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
    save_adata_outputAdded(data_version='split1', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
    save_adata_outputAdded(data_version='split2', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
    
