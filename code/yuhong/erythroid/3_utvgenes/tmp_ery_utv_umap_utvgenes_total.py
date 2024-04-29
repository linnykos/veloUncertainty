import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import bbknn

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
adata = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utv_preprocess.h5ad')

adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/erythroid_seed317_total_seurat.h5ad')
adata_total.var['highly_variable'] = adata.var['highly_variable'].copy()
adata_total.__dict__['_raw'].__dict__['_var'] = adata_total.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
total_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup6_eryutv_utvgenes/utv_erythroid_total_seurat_seed317.h5ad"
adata_total.write_h5ad(filename=total_path)
adata_total_res = utv.run_model(total_path, label, config_file=velo_config)
adata_total_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utvgenes_total.h5ad')
# using original adata umap
adata_total_res2 = adata_total_res.copy()
adata_total_res2.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(adata_total_res2,basis="umap",color=adata_total_res.uns['label'],
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_preumap_total.png")
# uncorrected recovered umap
sc.tl.umap(adata_total_res)
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color=label,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_uncorrected_total.png")
### batch correction
bbknn.bbknn(adata_total_res, batch_key='sequencing.batch')
adata_total_res.X = adata_total_res.X.toarray()
bbknn.ridge_regression(adata_total_res, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_total_res)
bbknn.bbknn(adata_total_res, batch_key='sequencing.batch')
sc.pp.neighbors(adata_total_res, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total_res)
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color=label,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_bbknn_total.png")

print("************* write adata_total_res to '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/' *******************")

