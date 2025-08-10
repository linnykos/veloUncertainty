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

def run_velovi_ery_Mark(data_version, method_prefix, gene_set_name, split_seed):
    split_version = data_version
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method = method_prefix + '_' + gene_set_name
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'
    from velovi import preprocess_data, VELOVI
    adata = None
    print_message_with_time("#################### "+dataset_long+'_'+method+'_'+data_version+": Read data ")
    adata = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_'+data_version+'.h5ad')
    print_message_with_time("#################### "+data_version+": Preprocess data ")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if (not 'wop' in method): adata = preprocess_data(adata)
    # train and apply model
    print_message_with_time("#################### "+data_version+": Train model ")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    # save vae
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    print_message_with_time("#################### "+data_version+": Save vae ")
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'.pt',overwrite=True)
    # write data
    print_message_with_time("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'.h5ad')
    print_message_with_time("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)

for i in range(4):
    split_seed = [320, 323, 326, 329][i]
    grid_seed = 227
    gene_set_name = 'nMark' + str(grid_seed)
    run_velovi_ery_Mark(data_version='total', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
    run_velovi_ery_Mark(data_version='split1', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
    run_velovi_ery_Mark(data_version='split2', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
