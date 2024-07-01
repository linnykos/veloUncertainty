import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

import matplotlib.pyplot as plt
import seaborn as sns

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_velovi/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/velovi/"
# load and preprocess data
adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad")
spliced = adata.layers['spliced'].copy() 
unspliced = adata.layers['unspliced'].copy()
gene_names = adata.var['Accession'].index.copy()

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# adata = preprocess_data(adata) 

# train and apply model
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()
vae.save(data_folder+'vae_erythroid_preprocess.pt')
print("########################### Wrote vae ###########################")

def add_velovi_outputs_to_adata(adata, vae):
    latent_time = vae.get_latent_time(n_samples=25)
    velocities = vae.get_velocity(n_samples=25, velo_statistic="mean")
    t = latent_time
    scaling = 20 / t.max(0)
    adata.layers["velocity"] = velocities / scaling
    adata.layers["latent_time_velovi"] = latent_time
    adata.var["fit_alpha"] = vae.get_rates()["alpha"] / scaling
    adata.var["fit_beta"] = vae.get_rates()["beta"] / scaling
    adata.var["fit_gamma"] = vae.get_rates()["gamma"] / scaling
    adata.var["fit_t_"] = (
        torch.nn.functional.softplus(vae.module.switch_time_unconstr)
        .detach()
        .cpu()
        .numpy()
    ) * scaling
    adata.layers["fit_t"] = latent_time.values * np.array(scaling)[np.newaxis, :] # scaling[np.newaxis, :] 
    adata.var['fit_scaling'] = 1.0

add_velovi_outputs_to_adata(adata, vae)
# ValueError: Multi-dimensional indexing (e.g. `obj[:, None]`) is no longer supported. 
### Convert to a numpy array before indexing instead.

scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', save=fig_folder+"erythroid_umap_preprocess.png")

# write data
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['highly_variable_genes'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# ValueError: '_index' is a reserved name for dataframe columns.
adata.write(filename=data_folder+"erythroid_velovi_preprocess.h5ad")
print("########################### Wrote data ###########################")

