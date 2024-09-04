split_seed = 317
dataset_long = 'larry'
dataset_short = 'larry'
method = 'utv'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_larry/'
data_version = 'split2'
celltype_label = 'state_info'

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_utv import *


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


###########
#total = sc.read_h5ad(data_folder+'v4_larry/larry_total_allgenes.h5ad')
#adata_split1 = sc.read_h5ad(data_folder+'v4_larry/seed317/seed317_larry_split1_allgenes.h5ad')
#adata_split2 = sc.read_h5ad(data_folder+'v4_larry/seed317/seed317_larry_split2_allgenes.h5ad')

utv_run_model_v4(data_version=data_version,method=method,dataset=dataset_short,velo_config=velo_config,data_folder=data_folder,split_seed=split_seed,celltype_label=celltype_label)


