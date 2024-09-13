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

def utv_larryMult(split_seed,data_version):
    dataset_long = 'larryMult'
    dataset_short = 'larryMult'
    method = 'utv'
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    celltype_label = 'state_info'

    data_path = data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_'+data_version+'_allgenes.h5ad'
    split = sc.read_h5ad(data_path)
    split.layers['spliced_original'] = split.layers['spliced'].copy()
    split.layers['unspliced_original'] = split.layers['unspliced'].copy()
    s = sc.read_h5ad(data_path)
    scv.pp.normalize_per_cell(s)
    scv.pp.log1p(s)
    sc.pp.highly_variable_genes(s)
    split.var['highly_variable'] = s.var['highly_variable'].copy()
    data_path_tmp = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/tmp/'+'v4_utv_'+dataset_short+'_'+data_version+'_'+str(split_seed)+'.h5ad'
    split.write(data_path_tmp)
    ### fit model
    split = utv.run_model(data_path_tmp, celltype_label, config_file=velo_config)
    print_message_with_time("#################### Write data ")
    split.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_utv_'+data_version+'_v4.h5ad')
    print_message_with_time("#################### All done for seed="+str(split_seed)+dataset_short+'+'+method+' '+data_version)

utv_larryMult(split_seed=323,data_version='split1')
utv_larryMult(split_seed=323,data_version='split2')
utv_larryMult(split_seed=326,data_version='split1')
utv_larryMult(split_seed=326,data_version='split2')
utv_larryMult(split_seed=329,data_version='split1')
utv_larryMult(split_seed=329,data_version='split2')

