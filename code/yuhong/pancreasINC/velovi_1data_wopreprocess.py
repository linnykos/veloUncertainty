import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import VELOVI
import datetime

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
savedata_folder = data_folder+"v2_pancreasINC/velovi/wopreprocess/"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

########################################
########### split1
data_version = "split1"
print_message_with_time("#################### "+data_version+": Read data ")
adata = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_'+data_version+'_allgenes.h5ad') 

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("#################### "+data_version+": Preprocess data ")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### "+data_version+": Train model ")
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()

# save vae
print_message_with_time("#################### "+data_version+": Save vae ")
vae.save(savedata_folder+'vae_panINC_velovi_'+data_version+'_wopreprocess_v2.pt',overwrite=True)

# save original counts
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
adata.layers['spliced_original'] = S_mat[:,positions]
adata.layers['unspliced_original'] = U_mat[:,positions]

# write data
print_message_with_time("#################### "+data_version+": Save adata (final version) ")
adata.write(filename=savedata_folder+"adata_panINC_velovi_"+data_version+"_wopreprocess_v2.h5ad")
print_message_with_time("#################### "+data_version+": All done for "+data_version)

########################################
########### split1
data_version = "split2"
print_message_with_time("#################### "+data_version+": Read data ")
adata = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_'+data_version+'_allgenes.h5ad') 

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("#################### "+data_version+": Preprocess data ")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### "+data_version+": Train model ")
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()
# save vae
print_message_with_time("#################### "+data_version+": Save vae ")
vae.save(savedata_folder+'vae_panINC_velovi_'+data_version+'_wopreprocess_v2.pt',overwrite=True)
# save original counts
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
adata.layers['spliced_original'] = S_mat[:,positions]
adata.layers['unspliced_original'] = U_mat[:,positions]
# write data
print_message_with_time("#################### "+data_version+": Save adata (final version) ")
adata.write(filename=savedata_folder+"adata_panINC_velovi_"+data_version+"_wopreprocess_v2.h5ad")
print_message_with_time("#################### "+data_version+": All done for "+data_version)

########################################
########### total
data_version = "total"
print_message_with_time("#################### "+data_version+": Read data ")
adata = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_'+data_version+'_allgenes.h5ad') 

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("#################### "+data_version+": Preprocess data ")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### "+data_version+": Train model ")
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()
# save vae
print_message_with_time("#################### "+data_version+": Save vae ")
vae.save(savedata_folder+'vae_panINC_velovi_'+data_version+'_wopreprocess_v2.pt',overwrite=True)
# save original counts
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
adata.layers['spliced_original'] = S_mat[:,positions]
adata.layers['unspliced_original'] = U_mat[:,positions]
# write data
print_message_with_time("#################### "+data_version+": Save adata (final version) ")
adata.write(filename=savedata_folder+"adata_panINC_velovi_"+data_version+"_wopreprocess_v2.h5ad")
print_message_with_time("#################### "+data_version+": All done for "+data_version)

