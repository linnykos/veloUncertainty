import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import datetime

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

def utv_run_model_v4(data_version,dataset,method,velo_config,data_folder,split_seed,celltype_label):
    data_path = None
    print_message_with_time("#################### Read data "+dataset+'+'+method+' '+data_version)
    if 'split' in data_version: 
        data_path = data_folder+'seed'+str(split_seed)+'_'+dataset+'_'+data_version+'_allgenes.h5ad'
    elif data_version == 'total':
        data_path = data_folder+dataset+'_total_allgenes.h5ad'
    adata = sc.read_h5ad(data_path)
    gene_names = adata.var.index.copy()
    S_mat = adata.layers['spliced'].copy()
    U_mat = adata.layers['unspliced'].copy()
    ### fit model
    print_message_with_time("#################### Fit model")
    adata = utv.run_model(data_path,celltype_label, config_file=velo_config)
    ### 
    positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
    positions = [positions_dict[gene] for gene in adata.var.index]
    adata.layers['spliced_original'] = S_mat[:,positions] # looks like i did not do this actually
    adata.layers['unspliced_original'] = U_mat[:,positions]
    print_message_with_time("#################### Write data ")
    adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset+'_utv_'+data_version+'_v4.h5ad')
    print_message_with_time("#################### All done for "+dataset+'+'+method+' '+data_version)
