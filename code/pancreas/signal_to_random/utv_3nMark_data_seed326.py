split_seed = 326
grid_seed = 236
gene_set_name = 'nMark'+str(grid_seed)
method = 'utv' + gene_set_name

dataset_long = 'pancreas'
dataset_short = 'pan'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
celltype_label = 'clusters'

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_utv import *

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

print('################ Write utv0 files #################')
def utv_run_model_nMark(data_version,velo_config,split_seed,grid_seed,celltype_label=celltype_label):
    gene_set_name = 'nMark'+str(grid_seed)
    method = 'utv_'+gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_pancreas/'
    adata = sc.read_h5ad(data_folder+dataset_short+'_'+data_version+'_'+gene_set_name+'.h5ad')
    adata.layers['spliced_original'] = adata.layers['spliced'].copy()
    adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
    adata.var['highly_variable'] = True
    adata.write_h5ad(data_folder+'utv0_'+dataset_short+'_'+data_version+'_'+gene_set_name+'.h5ad')
    adata = utv.run_model(data_folder+'utv0_'+dataset_short+'_'+data_version+'_'+gene_set_name+'.h5ad',
                            celltype_label, config_file=velo_config)
    adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+data_version+'_v4.h5ad')

data_version = 'split1'
utv_run_model_nMark(data_version, velo_config, split_seed, grid_seed,celltype_label=celltype_label)
data_version = 'split2'
utv_run_model_nMark(data_version, velo_config, split_seed, grid_seed,celltype_label=celltype_label)
data_version = 'total'
utv_run_model_nMark(data_version, velo_config, split_seed, grid_seed,celltype_label=celltype_label)
