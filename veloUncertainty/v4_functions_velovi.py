import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import datetime
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')


def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

def velovi_run_model(data_version,dataset,method,data_folder,split_seed):
    from velovi import preprocess_data, VELOVI
    adata = None
    print_message_with_time("#################### "+dataset+'_'+method+'_'+data_version+": Read data ")
    if 'split' in data_version: 
        adata = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_larry_'+data_version+'_allgenes.h5ad')
    elif data_version == 'total':
        adata = sc.read_h5ad(data_folder+dataset+'_total_allgenes.h5ad')
    gene_names = adata.var.index.copy()
    S_mat = adata.layers['spliced'].copy()
    U_mat = adata.layers['unspliced'].copy()
    print_message_with_time("#################### "+data_version+": Preprocess data ")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if 'wop' in method: adata = preprocess_data(adata)
    # train and apply model
    print_message_with_time("#################### "+data_version+": Train model ")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    # save vae
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    print_message_with_time("#################### "+data_version+": Save vae ")
    vae.save(savedata_folder+'vae_'+dataset+'_'+method+'_'+data_version+'_v4.pt',overwrite=True)
    # save original counts
    positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
    positions = [positions_dict[gene] for gene in adata.var.index]
    adata.layers['spliced_original'] = S_mat[:,positions]
    adata.layers['unspliced_original'] = U_mat[:,positions]
    # write data
    print_message_with_time("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset+'_'+method+'_'+data_version+'_v4.h5ad')
    print_message_with_time("#################### "+data_version+": All done for "+dataset+'_'+method+'_'+data_version)
