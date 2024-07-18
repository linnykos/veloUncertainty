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

# see:
## tmp_ery_utv_umap_utvgenes_total.py
## tmp_ery_utv_umap_utvgenes_317s1.py
## tmp_ery_utv_umap_utvgenes_317s2.py
## tmp_ery_utv_umap_utvgenes_320s1.py
## tmp_ery_utv_umap_utvgenes_320s2.py
