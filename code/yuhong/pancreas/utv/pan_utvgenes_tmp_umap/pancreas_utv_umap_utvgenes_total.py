import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
from sklearn.metrics.pairwise import cosine_similarity

### did not run pp.neighbors
# the below script uses the environment: "utvClone"

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = False
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.ASSIGN_POS_U = True

os.environ["TF_USE_LEGACY_KERAS"]="1"
# https://github.com/tensorflow/tensorflow/issues/62337
# To solve the following error: 
## ImportError: `keras.optimizers.legacy` is not supported in Keras 3. 
## When using `tf.keras`, to continue using a `tf.keras.optimizers.legacy` optimizer, 
## you can install the `tf_keras` package (Keras 2) and set the environment variable `TF_USE_LEGACY_KERAS=True` 
## to configure TensorFlow to use `tf_keras` when accessing `tf.keras`.
# Tried this but did not work:
## velo_config.TF_USE_LEGACY_KERAS=True

label='clusters'
adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utv_preprocess.h5ad")

adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_total_seurat.h5ad')
adata_total.var['highly_variable'] = adata.var['highly_variable'].copy()
adata_total.__dict__['_raw'].__dict__['_var'] = adata_total.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# To solve the error: ValueError: '_index' is a reserved name for dataframe columns.
# ref: https://github.com/theislab/scvelo/issues/255
total_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup5_panutv_utvgenes/tmp_pan_utvgenes_total_seurat_seed317.h5ad"
adata_total.write_h5ad(filename=total_path)
adata_total_res = utv.run_model(total_path, label, config_file=velo_config)
sc.tl.umap(adata_total_res)
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/utvgenes_total.png")
print("total counts done!")

adata_total_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_total.h5ad')

