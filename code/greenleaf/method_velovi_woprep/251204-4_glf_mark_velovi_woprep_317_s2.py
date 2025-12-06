split_seed = 317
grid_seed = 227

import scanpy as sc
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

#total_mark = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/glf_total_mark_genes.h5ad')
#split1_mark = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split1_mark_genes.h5ad')
split2_mark = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split2_mark_genes.h5ad')

method = 'velovi_woprep'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_velovi import *

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'

def velovi_run_model_glf_mark(adata, data_version,dataset_long,dataset_short,method,data_folder,split_seed):
    from velovi import preprocess_data, VELOVI
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
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/glf_mark/'
    print_message_with_time("#################### "+data_version+": Save vae ")
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'_glf_mark.pt',overwrite=True)
    # write data
    print_message_with_time("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_glf_mark.h5ad')
    print_message_with_time("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)


velovi_run_model_glf_mark(adata=split2_mark, data_version='split2', dataset_short=dataset_short, dataset_long=dataset_long, method=method,
						  data_folder=data_folder, split_seed=split_seed)


