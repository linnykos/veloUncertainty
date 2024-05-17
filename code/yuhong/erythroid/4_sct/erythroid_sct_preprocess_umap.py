import sctour as sct
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_sct/ery_sct_preprocess.h5ad")

# Visualization
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10, 9))
sc.pl.umap(adata, color='celltype', ax=axs[0, 0], legend_loc='on data', show=False, frameon=False)
sc.pl.umap(adata, color='sequencing.batch', ax=axs[0, 1], show=False, frameon=False)
sc.pl.umap(adata, color='ptime', ax=axs[1, 0], show=False, frameon=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', 
                         show=False, ax=axs[1, 1], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                         save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/sctour/erythroid_sct_vf_umappre.png")
# using original UMAP

## generate a UMAP embedding based on the inferred latent space
## may order the cells according to their pseudotime before this step, may yield a better trajectory for some datasets
adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['celltype', 'sequencing.batch'], legend_loc='on data', save="sct_umap2.png")


## largely unaffected by batch effects
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10, 9))
sc.pl.umap(adata, color='celltype', ax=axs[0, 0], legend_loc='on data', show=False, frameon=False)
sc.pl.umap(adata, color='sequencing.batch', ax=axs[0, 1], show=False, frameon=False)
sc.pl.umap(adata, color='ptime', ax=axs[1, 0], show=False, frameon=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', 
                         show=False, ax=axs[1, 1], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                         save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/sctour/erythroid_sct_vf_umap2000.png")
# Using newly computed UMAP (after sct training)











