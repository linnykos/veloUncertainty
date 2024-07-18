import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

import matplotlib.pyplot as plt
import seaborn as sns

# load and preprocess data
adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad")
spliced = adata.layers['spliced'].copy()
unspliced = adata.layers['unspliced'].copy()
gene_names = adata.var['highly_variable_genes'].copy()

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# adata = preprocess_data(adata) # 3696 × 1074

# train and apply model
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()

# save vae
vae.save('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_velovi/vae_preprocess.pt',overwrite=True)
print("########################### vae saved ###########################")

# write adata
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['highly_variable_genes'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset
adata.write(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_velovi/pancreas_velovi_preprocess.h5ad")
print("########################### adata saved ###########################")
