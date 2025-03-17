import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_utv import *
from v4_functions import *


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




#fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'+'seed'+str(split_seed)+'/'+method+'/'

data_version = 'split1'
split_seed = 317

def utv_run_model_ery_Mark(gene_set_name, split_seed, data_version, velo_config):
    dataset_long = 'erythroid'
    dataset_short = 'ery'
    method_prefix = 'utv'
    method = method_prefix + '_' + gene_set_name
    celltype_label = get_celltype_label(dataset_short)
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    data_path = None
    print_message_with_time("#################### Read data "+dataset_short+'+'+method+' '+data_version)
    if 'split' in data_version: 
        data_path = data_folder+'utv0_adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_'+data_version+'_allgenes.h5ad'
        #data_folder+'seed'+str(split_seed)+'_'+dataset+'_'+data_version+'_allgenes.h5ad'
    elif data_version == 'total':
        data_path = data_folder+'utv0_adata_'+dataset_short+'_'+gene_set_name+'_total_allgenes.h5ad'
        #data_path = data_folder+dataset+'_total_allgenes.h5ad'
    adata = sc.read_h5ad(data_path)
    gene_names = adata.var.index.copy()
    ### fit model
    print_message_with_time("#################### Fit model")
    adata = utv.run_model(data_path, celltype_label, config_file=velo_config)
    ### 
    print_message_with_time("#################### Write data ")
    adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+data_version+'.h5ad')
    print_message_with_time("#################### All done for "+dataset_short+'+'+method+' '+data_version)

utv_run_model_ery_Mark(gene_set_name='Mark', split_seed=317, data_version='total', velo_config=velo_config)
utv_run_model_ery_Mark(gene_set_name='Mark', split_seed=317, data_version='split1', velo_config=velo_config)
utv_run_model_ery_Mark(gene_set_name='Mark', split_seed=317, data_version='split2', velo_config=velo_config)
