split_seed = 320
grid_seed=230

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_velovi import *

dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

method = 'velovi_woprep'
file_suffix = 'glf_gseactrl'+str(grid_seed)

gene_set_name = 'gseactrl'

total = sc.read(data_folder+'v4_greenleaf/glf_genes_gsea/glf_total_'+gene_set_name+str(grid_seed)+'.h5ad')
split1 = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_gsea/glf_split1_'+gene_set_name+str(grid_seed)+'.h5ad')
split2 = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_gsea/glf_split2_'+gene_set_name+str(grid_seed)+'.h5ad')

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'

# data_suffix='gsea', 'gseactrl'+str(grid_seed)
def velovi_run_model_glf_gsea(adata, data_version,dataset_long,dataset_short,method,data_folder,split_seed, data_suffix):
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
    print_message_with_time("#################### "+data_version+": Finished training ")
    # save vae
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'_'+data_suffix+'/'
    print_message_with_time("#################### "+data_version+": Save vae ")
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'_'+file_suffix+'.pt',overwrite=True)
    # write data
    print_message_with_time("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_'+file_suffix+'.h5ad')
    print_message_with_time("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)
    print('#################### Data saved at: '+savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_'+file_suffix+'.h5ad')


velovi_run_model_glf_gsea(adata=total, data_version='total', dataset_short=dataset_short, dataset_long=dataset_long, method=method,
						  data_folder=data_folder, split_seed=split_seed, data_suffix='gseactrl'+str(grid_seed))
velovi_run_model_glf_gsea(adata=split1, data_version='split1', dataset_short=dataset_short, dataset_long=dataset_long, method=method,
						  data_folder=data_folder, split_seed=split_seed, data_suffix='gseactrl'+str(grid_seed))
velovi_run_model_glf_gsea(adata=split2, data_version='split2', dataset_short=dataset_short, dataset_long=dataset_long, method=method,
						  data_folder=data_folder, split_seed=split_seed, data_suffix='gseactrl'+str(grid_seed))



