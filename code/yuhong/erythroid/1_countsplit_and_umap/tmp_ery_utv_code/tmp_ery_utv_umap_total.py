import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os

# the below script uses the environment: "utvClone"

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1

os.environ["TF_USE_LEGACY_KERAS"]="1"
# https://github.com/tensorflow/tensorflow/issues/62337
# To solve the following error: 
## ImportError: `keras.optimizers.legacy` is not supported in Keras 3. 
## When using `tf.keras`, to continue using a `tf.keras.optimizers.legacy` optimizer, 
## you can install the `tf_keras` package (Keras 2) and set the environment variable `TF_USE_LEGACY_KERAS=True` 
## to configure TensorFlow to use `tf_keras` when accessing `tf.keras`.
# Tried this but did not work:
## velo_config.TF_USE_LEGACY_KERAS=True

label='celltype'
adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad")
#adata = scv.datasets.gastrulation_erythroid()
scv.pp.filter_genes(adata, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) 
### Extracted 2000 highly variable genes.
scv.pp.log1p(adata)

true_indices = adata.var['highly_variable'][adata.var['highly_variable'] == True].index.tolist()
velo_config.VGENES = true_indices

adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_total_seurat.h5ad')
adata_total.var['highly_variable'] = adata.var['highly_variable'].copy()
adata_total.__dict__['_raw'].__dict__['_var'] = adata_total.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# To solve the error: ValueError: '_index' is a reserved name for dataframe columns.
# ref: https://github.com/theislab/scvelo/issues/255
total_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup3_eryutv/out/utv_erythroid_total_seurat_seed317.h5ad"
adata_total.write_h5ad(filename=total_path)
adata_total_res = utv.run_model(total_path, label, config_file=velo_config)
adata_total_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_total.h5ad')
sc.tl.umap(adata_total_res)
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color=label,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_uncorrected_total.png")
# using original adata umap
adata_total_res.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color=adata_total_res.uns['label'],
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_preumap_total.png")
print("************* write adata_total_res to '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_total.h5ad' *******************")

