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
adata_total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")

gene_names = adata_total.var.index.copy()
S_mat_total = adata_total.layers['spliced'].copy()
U_mat_total = adata_total.layers['unspliced'].copy()

print_message_with_time("#################### Preprocess data ")
scv.pp.filter_and_normalize(adata_total, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
adata_total = preprocess_data(adata_total)

# train and apply model
print_message_with_time("#################### Train model ")
VELOVI.setup_anndata(adata_total, spliced_layer="Ms", unspliced_layer="Mu")
vae_total = VELOVI(adata_total)
vae_total.train()

# save vae
print_message_with_time("#################### Save vae ")
vae_total.save(data_folder+'v2_erythroid/velovi/vae_ery_velovi_total_v2.pt',overwrite=True)
adata_total.write(filename=data_folder+"v2_erythroid/velovi/backup/adata_ery_velovi_total.h5ad")

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata_total.var.index]

S_mat_total = S_mat_total[:,positions]
U_mat_total = U_mat_total[:,positions]
adata_total.layers['spliced_original'] = S_mat_total
adata_total.layers['unspliced_original'] = U_mat_total

# write data
print_message_with_time("#################### Save adata (final version) ")
adata_total.write(filename=data_folder+"v2_erythroid/velovi/adata_ery_velovi_total_v2.h5ad")
print_message_with_time("#################### All done for total ")
