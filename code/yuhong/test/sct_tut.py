import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/EX/EX_development_human_cortex_10X.h5ad")
adata.shape  # (36318,19073)

#adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad")

sc.pl.umap(adata, color=['celltype', 'Sample batch'], legend_loc='on data', save="sct_umap1.png")

# Model training
## count the number of genes detected in each cell
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
## select highly variable genes
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=1000, subset=True)  # 36318 × 1000
# ImportError: Please install skmisc package via `pip install --user scikit-misc

# default loss_model uses NB conditioned likelihood on raw UMI counts (in adata.X)
# default "percent" of cells used to train the model is set to 0.9 when #total_cells<10000, =0.2 when #total_cells>10000
# alpha_recon_lec, alpha_recon_lode: default 0.5 (balance reconstruction errors from encoder-derived latent space and ODE-solver-derived latent space)
tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
tnode.train() # 14:20 ~ 14:49

# Infer cellular dynamics
## pseudotime
adata.obs['ptime'] = tnode.get_time()

## latent space
### alpha_z and alpha_predz adjust weights given to the latent space from variational infereence and that from ODE solver
### larger alpha_z skews the latent space towards the intrinsic transcriptomic structure
### larger alpha_predz is more representative of the extrinsic pseudotime ordering
#zs represents the latent z from variational inference, and pred_zs represents the latent z from ODE solver
#mix_zs represents the weighted combination of the two, which is used for downstream analysis
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs
## adata.obsm['X_TNODE'].shape = (36318, 5)

## vector field
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
## adata.obsm['X_VF'].shape = (36318, 5)

# Visualization
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))
sc.pl.umap(adata, color='celltype', ax=axs[0, 0], legend_loc='on data', show=False, frameon=False)
sc.pl.umap(adata, color='Sample batch', ax=axs[0, 1], show=False, frameon=False)
sc.pl.umap(adata, color='ptime', ax=axs[1, 0], show=False, frameon=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', 
                         color='celltype', show=False, ax=axs[1, 1], legend_loc='none', frameon=False, 
                         size=100, alpha=0.2, save="../../../git/veloUncertainty/fig/yuhong/sct_tut/sct_vf_umap1.png")
# using original UMAP

## generate a UMAP embedding based on the inferred latent space
## may order the cells according to their pseudotime before this step, may yield a better trajectory for some datasets
adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=15)
sc.tl.umap(adata, min_dist=0.1)
sc.pl.umap(adata, color=['celltype', 'Sample batch'], legend_loc='on data', save="sct_umap2.png")


## largely unaffected by batch effects
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))
sc.pl.umap(adata, color='celltype', ax=axs[0, 0], legend_loc='on data', show=False, frameon=False)
sc.pl.umap(adata, color='Sample batch', ax=axs[0, 1], show=False, frameon=False)
sc.pl.umap(adata, color='ptime', ax=axs[1, 0], show=False, frameon=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', 
                         color='celltype', show=False, ax=axs[1, 1], legend_loc='none', frameon=False, 
                         size=100, alpha=0.2, save="../../../git/veloUncertainty/fig/yuhong/sct_tut/sct_vf_umap2.png")
# Using newly computed UMAP (after sct training)











