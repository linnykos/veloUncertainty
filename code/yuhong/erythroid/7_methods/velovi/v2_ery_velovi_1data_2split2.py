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

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

print_message_with_time("#################### Read data ")
adata_split2 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split2_allgenes.h5ad')

gene_names = adata_split2.var.index.copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()

print_message_with_time("#################### Preprocess data ")
scv.pp.filter_and_normalize(adata_split2, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)
adata_split2 = preprocess_data(adata_split2)

# train and apply model
print_message_with_time("#################### Train model ")
VELOVI.setup_anndata(adata_split2, spliced_layer="Ms", unspliced_layer="Mu")
vae_split2 = VELOVI(adata_split2)
vae_split2.train()

# save vae
print_message_with_time("#################### Save vae ")
vae_split2.save(data_folder+'v2_erythroid/velovi/vae_ery_velovi_split2_v2.pt',overwrite=True)
adata_split2.write(filename=data_folder+"v2_erythroid/velovi/backup/adata_ery_velovi_split2.h5ad")

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata_split2.var.index]

S_mat_split2 = S_mat_split2[:,positions]
U_mat_split2 = U_mat_split2[:,positions]
adata_split2.layers['spliced_original'] = S_mat_split2
adata_split2.layers['unspliced_original'] = U_mat_split2

# write data
print_message_with_time("#################### Save adata (final version) ")
adata_split2.write(filename=data_folder+"v2_erythroid/velovi/adata_ery_velovi_split2_v2.h5ad")
print_message_with_time("#################### All done for split2 ")
