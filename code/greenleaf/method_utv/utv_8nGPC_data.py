split_seed = 317
dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'utv_nGPC227'
#data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'
celltype_label = 'cluster_name'

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os

# the below script uses the environment: "utvClone"
# larry/erythroid configuration
velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True # default value
# (bool) linear regression $R^2$ on extreme quantile (default) or full data (adjusted)
# valid when self.VGENES = 'basic'
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1  # default initialization
# (float) threshold of R2 at later stage of the optimization proces
# to capture the dynamics of more genes beside initially selected velocity genes
# self.AGENES_R2 = 1 will switch to origianl mode with no amplification
os.environ["TF_USE_LEGACY_KERAS"]="1"

print('################ Write utv0 files #################')
def utv_run_model_v4_nGPC227(data_version,velo_config,data_folder,split_seed,method,celltype_label='cluster_name'):
    if data_version=='total':
        adata = sc.read_h5ad(data_folder+'glf_'+data_version+'_'+method+'_utv0.h5ad')
        adata.layers['spliced_original'] = adata.layers['spliced'].copy()
        adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
        adata.var['highly_variable'] = True
        adata.write_h5ad(data_folder+'glf_'+data_version+'_'+method+'_utv0.h5ad')
        adata = utv.run_model(data_folder+'glf_'+data_version+'_'+method+'_utv0.h5ad',celltype_label, config_file=velo_config)
        adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_glf_'+method+'_'+data_version+'_v4.h5ad')
    elif 'split' in data_version:
        adata = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_glf_'+data_version+'_'+method+'_utv0.h5ad')
        adata.layers['spliced_original'] = adata.layers['spliced'].copy()
        adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
        adata.var['highly_variable'] = True
        adata.write_h5ad(data_folder+'seed'+str(split_seed)+'_glf_'+data_version+'_'+method+'_utv0.h5ad')
        adata = utv.run_model(data_folder+'seed'+str(split_seed)+'_glf_'+data_version+'_'+method+'_utv0.h5ad',
                              celltype_label, config_file=velo_config)
        adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_glf_'+method+'_'+data_version+'_v4.h5ad')

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'

# total
utv_run_model_v4_nGPC227(data_version='total', velo_config=velo_config, data_folder=data_folder, split_seed=317, method=method,celltype_label='cluster_name')
# split1
#utv_run_model_v4_nGPC227(data_version='split1', velo_config=velo_config, data_folder=data_folder, split_seed=317, method=method,celltype_label='cluster_name')
# split2
#utv_run_model_v4_nGPC227(data_version='split2', velo_config=velo_config, data_folder=data_folder, split_seed=317, method=method,celltype_label='cluster_name')

