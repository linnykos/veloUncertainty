import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
savedata_folder = data_folder+"v2_erythroid/velovi/preprocessed/"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

############################################
### split1
print_message_with_time("#################### split1: Read data ")
adata_split1 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split1_allgenes.h5ad')

gene_names = adata_split1.var.index.copy()
S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()

print_message_with_time("#################### split1: Preprocess data ")
scv.pp.filter_and_normalize(adata_split1, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata_split1, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### split1: Train model ")
VELOVI.setup_anndata(adata_split1, spliced_layer="Ms", unspliced_layer="Mu")
vae_split1 = VELOVI(adata_split1)
vae_split1.train()

# save vae
print_message_with_time("#################### split1: Save vae ")
vae_split1.save(savedata_folder+'vae_ery_velovi_split1_preprocessed_v2.pt',overwrite=True)

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata_split1.var.index]

S_mat_split1 = S_mat_split1[:,positions]
U_mat_split1 = U_mat_split1[:,positions]
adata_split1.layers['spliced_original'] = S_mat_split1
adata_split1.layers['unspliced_original'] = U_mat_split1

# write data
print_message_with_time("#################### split1: Save adata (final version) ")
adata_split1.write(filename=savedata_folder+"adata_ery_velovi_split1_preprocessed_v2.h5ad")
print_message_with_time("#################### split1: All done for split1 ")

############################################
### split2
print_message_with_time("#################### split2: Read data ")
adata_split2 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split2_allgenes.h5ad')

gene_names = adata_split2.var.index.copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()

print_message_with_time("#################### split2: Preprocess data ")
scv.pp.filter_and_normalize(adata_split2, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### split2: Train model ")
VELOVI.setup_anndata(adata_split2, spliced_layer="Ms", unspliced_layer="Mu")
vae_split2 = VELOVI(adata_split2)
vae_split2.train()

# save vae
print_message_with_time("#################### split2: Save vae ")
vae_split2.save(savedata_folder+'vae_ery_velovi_split2_preprocessed_v2.pt',overwrite=True)

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata_split2.var.index]

S_mat_split2 = S_mat_split2[:,positions]
U_mat_split2 = U_mat_split2[:,positions]
adata_split2.layers['spliced_original'] = S_mat_split2
adata_split2.layers['unspliced_original'] = U_mat_split2

# write data
print_message_with_time("#################### Save adata (final version) ")
adata_split2.write(filename=savedata_folder+"adata_ery_velovi_split2_preprocessed_v2.h5ad")
print_message_with_time("#################### All done for split2 ")

############################################
### total
print_message_with_time("#################### total: Read data ")
adata_total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")

gene_names = adata_total.var.index.copy()
S_mat_total = adata_total.layers['spliced'].copy()
U_mat_total = adata_total.layers['unspliced'].copy()

print_message_with_time("#################### total: Preprocess data ")
scv.pp.filter_and_normalize(adata_total, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### total: Train model ")
VELOVI.setup_anndata(adata_total, spliced_layer="Ms", unspliced_layer="Mu")
vae_total = VELOVI(adata_total)
vae_total.train()

# save vae
print_message_with_time("#################### total: Save vae ")
vae_total.save(savedata_folder+'vae_ery_velovi_total_preprocessed_v2.pt',overwrite=True)

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata_total.var.index]

S_mat_total = S_mat_total[:,positions]
U_mat_total = U_mat_total[:,positions]
adata_total.layers['spliced_original'] = S_mat_total
adata_total.layers['unspliced_original'] = U_mat_total

# write data
print_message_with_time("#################### total: Save adata (final version) ")
adata_total.write(filename=savedata_folder+"adata_ery_velovi_total_preprocessed_v2.h5ad")
print_message_with_time("#################### total: All done for total ")

