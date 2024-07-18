import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

import matplotlib.pyplot as plt
import seaborn as sns

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_velovi/"
s1_317_name = "pancreas_seed317_split1_seurat.h5ad"
s2_317_name = "pancreas_seed317_split2_seurat.h5ad"
total_name = "pancreas_velovi_preprocess.h5ad"

figure_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/velovi/"

total = sc.read(data_folder+total_name)

label = 'clusters'
label_color = 'clusters_colors'

## scatter plot without velocity arrows
scv.pl.scatter(total, color=label, cmap=label_color, save=figure_folder+"pancreas_preprocess_scatter.png")

def train_vae_and_save(adata,adata_name,vae_name):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    adata.write(filename=data_folder+adata_name+".h5ad")
    print("########################### adata splits saved ###########################")
    vae.save(data_folder+vae_name+'.pt',overwrite=True)



# 317 split1
s1_317 = sc.read(data_folder+s1_317_name)
train_vae_and_save(adata=s1_317, adata_name="pancreas_seed317_split1", vae_name="vae_seed317_split1")

# 317 split2
s2_317 = sc.read(data_folder+s2_317_name)
train_vae_and_save(adata=s2_317, adata_name="pancreas_seed317_split2", vae_name="vae_seed317_split2")

#s1_317.write(filename=data_folder+"pancreas_test.h5ad")
#vae_317s1.save(data_folder+'vae_test.pt',overwrite=True)
#test_adata = sc.read(data_folder+"pancreas_test.h5ad")
#test_vae = VELOVI.load(data_folder+'vae_test.pt',test_adata)





