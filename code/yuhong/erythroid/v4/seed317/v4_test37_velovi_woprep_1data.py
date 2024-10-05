
dataset_short = 'ery'
dataset_long = 'erythroid'
method = 'velovi_woprep'
split_seed = 317
celltype_label = 'celltype'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_test/v4_"+dataset_long+'/'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI


def velovi_run_model_v4test(data_version,dataset_long,dataset_short,method,data_folder,split_seed):
    from velovi import preprocess_data, VELOVI
    adata = None
    print("#################### "+dataset_long+'_'+method+'_'+data_version+": Read data ")
    adata = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_allgenes_'+data_version+'.h5ad')
    adata.layers['spliced_original'] = adata.layers['spliced'].copy()
    adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
    print("#################### "+data_version+": Preprocess data ")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if (not 'wop' in method): adata = preprocess_data(adata)
    # train and apply model
    print("#################### "+data_version+": Train model ")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    # save vae
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    print("#################### "+data_version+": Save vae ")
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'_v4.pt',overwrite=True)
    # write data
    print("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_v4.h5ad')
    print("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)

velovi_run_model_v4test(data_version='37-3',dataset_short=dataset_short,dataset_long=dataset_long,method=method,data_folder=data_folder,split_seed=split_seed)

velovi_run_model_v4test(data_version='37-7',dataset_short=dataset_short,dataset_long=dataset_long,method=method,data_folder=data_folder,split_seed=split_seed)
