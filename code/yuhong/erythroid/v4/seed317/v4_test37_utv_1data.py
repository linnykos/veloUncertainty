split_seed = 317
dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'utv'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_test/v4_'+dataset_long+'/'
data_version = 'ep37-3'
celltype_label = 'celltype'

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os

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
def utv_run_model_v4_test(data_version,dataset,method,velo_config,data_folder,split_seed,celltype_label):
    data_path = None
    print("#################### Read data "+dataset+'+'+method+' '+data_version)
    data_path = data_folder+'seed'+str(split_seed)+'_'+dataset+'_allgenes_'+data_version+'.h5ad'
    adata = sc.read_h5ad(data_path)
    adata.layers['spliced_original'] = adata.layers['spliced'].copy()
    adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
    ### fit model
    print("#################### Fit model")
    adata = utv.run_model(data_path,celltype_label, config_file=velo_config)
    ### 
    print("#################### Write data ")
    adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset+'_utv_'+data_version+'_v4.h5ad')
    print("#################### All done for "+dataset+'+'+method+' '+data_version)

utv_run_model_v4_test(data_version='37-3',method=method,dataset=dataset_short,velo_config=velo_config,data_folder=data_folder,split_seed=split_seed,celltype_label=celltype_label)
utv_run_model_v4_test(data_version='37-7',method=method,dataset=dataset_short,velo_config=velo_config,data_folder=data_folder,split_seed=split_seed,celltype_label=celltype_label)


