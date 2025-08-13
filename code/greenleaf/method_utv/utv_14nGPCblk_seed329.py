split_seed = 329
grid_seed = 239

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os

dataset_long = 'greenleaf'
dataset_short = 'glf'
celltype_label = 'cluster_name'

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
def utv_run_model_v4_nGPCblk(data_version,velo_config,split_seed,grid_seed,celltype_label='cluster_name'):
    gene_set_name = 'nGPCblk'+str(grid_seed)
    method = 'utv_'+gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'
    if data_version=='total':
        adata = sc.read_h5ad(data_folder+dataset_short+'_'+data_version+'_'+gene_set_name+'.h5ad')
        adata.layers['spliced_original'] = adata.layers['spliced'].copy()
        adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
        adata.var['highly_variable'] = True
        adata.write_h5ad(data_folder+'utv0_'+dataset_short+'_'+data_version+'_'+method+'.h5ad')
        adata = utv.run_model(data_folder+'utv0_'+dataset_short+'_'+data_version+'_'+method+'.h5ad',
                              celltype_label, config_file=velo_config)
        adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+data_version+'_v4.h5ad')
    elif 'split' in data_version:
        adata = sc.read_h5ad(data_folder+dataset_short+'_seed'+str(split_seed)+'_'+data_version+'_'+gene_set_name+'.h5ad')
        adata.layers['spliced_original'] = adata.layers['spliced'].copy()
        adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
        adata.var['highly_variable'] = True
        adata.write_h5ad(data_folder+'utv0_'+dataset_short+'_seed'+str(split_seed)+'_'+data_version+'_'+method+'.h5ad')
        adata = utv.run_model(data_folder+'utv0_'+dataset_short+'_seed'+str(split_seed)+'_'+data_version+'_'+method+'.h5ad',
                              celltype_label, config_file=velo_config)
        adata.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+data_version+'_v4.h5ad')


# total
utv_run_model_v4_nGPCblk(data_version='total', velo_config=velo_config, split_seed=split_seed, grid_seed=grid_seed, celltype_label='cluster_name')
# split1
utv_run_model_v4_nGPCblk(data_version='split1', velo_config=velo_config, split_seed=split_seed, grid_seed=grid_seed, celltype_label='cluster_name')
# split2
utv_run_model_v4_nGPCblk(data_version='split2', velo_config=velo_config, split_seed=split_seed, grid_seed=grid_seed, celltype_label='cluster_name')

