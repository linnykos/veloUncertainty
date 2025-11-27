import scanpy as sc

dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'velovi_woprep'
split_seed = 317
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

adata1 = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes_shrunk_disp.h5ad')
#adata2 = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes_shrunk_disp.h5ad')
#adata = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/'+dataset_short+'_total_allgenes_shrunk_disp.h5ad')

## run velovi
data_version = "split1"
celltype_label = 'celltype'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/shrunk_disp/'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_velovi import *


def velovi_run_model_shrunk_disp(adata, data_version,dataset_long,dataset_short,method,data_folder,split_seed):
    from velovi import preprocess_data, VELOVI
    adata = None
    print_message_with_time("#################### "+dataset_long+'_'+method+'_'+data_version+": Read data ")
    adata = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version,allgenes=True,outputAdded=False)
    gene_names = adata.var.index.copy()
    adata.layers['spliced_original'] = adata.layers['spliced'].copy()
    adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
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
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'_shrunk_disp.pt',overwrite=True)
    # write data
    print_message_with_time("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_shrunk_disp.h5ad')
    print_message_with_time("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)


velovi_run_model_shrunk_disp(adata=adata1, data_version=data_version,dataset_short=dataset_short,dataset_long=dataset_long,method=method,
	data_folder=data_folder,split_seed=split_seed)

