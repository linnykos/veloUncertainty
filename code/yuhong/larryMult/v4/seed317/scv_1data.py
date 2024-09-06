import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import read_data_v4
from v4_functions_scv import scv_compute_velocity

split_seed = 317
dataset_long = "larryMult"
dataset_short = "larryMult"
method = "scv"

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

print_message_with_time("################## Read data")
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=True,outputAdded=False)
adata_split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=True,outputAdded=False)
adata_split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=True,outputAdded=False)

adata_split1.layers['spliced_original'] = adata_split1.layers['spliced'].copy()
adata_split1.layers['unspliced_original'] = adata_split1.layers['unspliced'].copy()
adata_split2.layers['spliced_original'] = adata_split2.layers['spliced'].copy()
adata_split2.layers['unspliced_original'] = adata_split2.layers['unspliced'].copy()
total.layers['spliced_original'] = total.layers['spliced'].copy()
total.layers['unspliced_original'] = total.layers['unspliced'].copy()


print_message_with_time("################## Run model on total")
scv_compute_velocity(total,dataset_short) 

print_message_with_time("################## Read model on split1")
scv_compute_velocity(adata_split1,dataset_short)

print_message_with_time("################## Read model on split2")
scv_compute_velocity(adata_split2,dataset_short)

### add umapOriginal
total.obsm['X_umapOriginal'] = total.obsm['X_umap'].copy()
total.obsm['X_umapOriginal'][:,0] = np.array(total.obs['SPRING-x'])
total.obsm['X_umapOriginal'][:,1] = np.array(total.obs['SPRING-y'])

adata_split1.obsm['X_umapOriginal'] = adata_split1.obsm['X_umap'].copy()
adata_split1.obsm['X_umapOriginal'][:,0] = np.array(adata_split1.obs['SPRING-x'])
adata_split1.obsm['X_umapOriginal'][:,1] = np.array(adata_split1.obs['SPRING-y'])

adata_split2.obsm['X_umapOriginal'] = adata_split2.obsm['X_umap'].copy()
adata_split2.obsm['X_umapOriginal'][:,0] = np.array(adata_split2.obs['SPRING-x'])
adata_split2.obsm['X_umapOriginal'][:,1] = np.array(adata_split2.obs['SPRING-y'])


# write data
print_message_with_time("################## Write data")
total.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4.h5ad')
adata_split1.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
adata_split2.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')

print_message_with_time("################## All done with the data")

