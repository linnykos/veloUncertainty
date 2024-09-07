split_seed = 317
dataset_long = 'larryMult'
dataset_short = 'larryMult'
method = 'utv'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
data_version = 'split1'
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

"""
# try to use a subset of genes selected from total
## total -> filter out 5k genes
adata = sc.read_h5ad(data_folder+dataset_short+'_total_allgenes.h5ad')
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=5000) # try to filter 5k genes out, got 2831
adata.var.index
## split
del data_path
split1_subset = split1[:, split1.var_names.isin(adata.var.index)]
# 5k genes
data_path_tmp = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/tmp/'+'v4_utv_'+dataset_short+'_'+data_version+'_5k.h5ad'
split1_subset.write(data_path_tmp)
# 2k genes
data_path_tmp = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/tmp/'+'v4_utv_'+dataset_short+'_'+data_version+'_2k.h5ad'
split1_subset.write(data_path_tmp)
"""
data_path = data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_'+data_version+'_allgenes.h5ad'
split1 = sc.read_h5ad(data_path)

split1.layers['spliced_original'] = split1.layers['spliced'].copy()
split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()

s1 = sc.read_h5ad(data_path)
scv.pp.normalize_per_cell(s1)
scv.pp.log1p(s1)
sc.pp.highly_variable_genes(s1)

split1.var['highly_variable'] = s1.var['highly_variable'].copy()

data_path_tmp = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/tmp/'+'v4_utv_'+dataset_short+'_'+data_version+'.h5ad'
split1.write(data_path_tmp)

### fit model
split1 = utv.run_model(data_path_tmp,celltype_label, config_file=velo_config)
# Extracted 343 highly variable genes.
# Computing moments for 1112 genes with n_neighbors: 30 and n_pcs: 30
### 
print_message_with_time("#################### Write data ")
split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_utv_'+data_version+'_v4.h5ad')
print_message_with_time("#################### All done for "+dataset_short+'+'+method+' '+data_version)

