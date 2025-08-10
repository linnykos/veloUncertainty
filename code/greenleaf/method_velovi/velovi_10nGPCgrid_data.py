dataset_short = 'glf'
dataset_long = 'greenleaf'
gene_set_name = 'nGPCgrid'
method = 'velovi_'+gene_set_name
split_seed = 317
celltype_label = 'cluster_name'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

def velovi_run_model(adata,data_version,dataset_long,dataset_short,method,data_folder,split_seed):
    from velovi import preprocess_data, VELOVI
    gene_names = adata.var.index.copy()
    S_mat = adata.layers['spliced'].copy()
    U_mat = adata.layers['unspliced'].copy()
    adata.layers['spliced_original'] = S_mat
    adata.layers['unspliced_original'] = U_mat
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
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'_GPC.pt',overwrite=True)
    # save original counts
    positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
    positions = [positions_dict[gene] for gene in adata.var.index]
    adata.layers['spliced_original'] = S_mat[:,positions]
    adata.layers['unspliced_original'] = U_mat[:,positions]
    # write data
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_GPC.h5ad')
    print("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)

total = sc.read(data_folder+dataset_short+'_total_'+gene_set_name+'.h5ad')
split1 = sc.read(data_folder+dataset_short+'_seed'+str(split_seed)+'_split1_'+gene_set_name+'.h5ad')
split2 = sc.read(data_folder+dataset_short+'_seed'+str(split_seed)+'_split2_'+gene_set_name+'.h5ad')

data_version = 'split1'
velovi_run_model(adata=split1,data_version=data_version,dataset_long=dataset_long,dataset_short=dataset_short,method=method,data_folder=data_folder,split_seed=split_seed)

data_version = 'split2'
velovi_run_model(adata=split2,data_version=data_version,dataset_long=dataset_long,dataset_short=dataset_short,method=method,data_folder=data_folder,split_seed=split_seed)

data_version = 'total'
velovi_run_model(adata=total,data_version=data_version,dataset_long=dataset_long,dataset_short=dataset_short,method=method,data_folder=data_folder,split_seed=split_seed)

