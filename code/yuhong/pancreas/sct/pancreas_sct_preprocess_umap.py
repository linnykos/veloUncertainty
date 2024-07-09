import sctour as sct
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_sct/pan_sct_preprocess.h5ad")

# at location: /home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour
fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
sc.pl.umap(adata, color='clusters', ax=axs[0], legend_loc='on data', show=False, frameon=False)
sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='clusters', 
                         show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                         save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/pancreas_sct_vf_preumap.png")

adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['clusters'], legend_loc='on data', 
           save="pancreas_sct_umap2000.png") # manually move out of the figures/ folder, "umappancreas_sct_umap2000.png"

## largely unaffected by batch effects
fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))
sc.pl.umap(adata, color='clusters', ax=axs[0], legend_loc='on data', show=False, frameon=False)
sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='clusters', 
                         show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                         save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/pancreas_sct_vf_umap2000.png")







