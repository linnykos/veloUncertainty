dataset_long = "pancreas"
dataset_short = "pan"
method = "scv"
split_seed = 317

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_scv import scv_compute_velocity

total = read_data_v4(dataset_long,dataset_short,method,split_seed=split_seed,data_version='total',allgenes=True,outputAdded=False)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed=split_seed,data_version='split1',allgenes=True,outputAdded=False)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed=split_seed,data_version='split2',allgenes=True,outputAdded=False)

split1.layers['spliced_original'] = split1.layers['spliced'].copy()
split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()
split2.layers['spliced_original'] = split1.layers['spliced'].copy()
split2.layers['unspliced_original'] = split1.layers['unspliced'].copy()
total.layers['spliced_original'] = total.layers['spliced'].copy()
total.layers['unspliced_original'] = total.layers['unspliced'].copy()

scv_compute_velocity(total,dataset_short) 
scv_compute_velocity(split1,dataset_short) 
scv_compute_velocity(split2,dataset_short) 

total.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_total_v4.h5ad')
split1.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
split2.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
