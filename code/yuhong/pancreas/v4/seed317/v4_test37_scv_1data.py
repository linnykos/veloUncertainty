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

split1 = sc.read_h5ad(data_folder+'v4_test/v4_pancreas/seed317_pan_allgenes_37-3.h5ad')
split2 = sc.read_h5ad(data_folder+'v4_test/v4_pancreas/seed317_pan_allgenes_37-7.h5ad')

split1.layers['spliced_original'] = split1.layers['spliced'].copy()
split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()
split2.layers['spliced_original'] = split1.layers['spliced'].copy()
split2.layers['unspliced_original'] = split1.layers['unspliced'].copy()

scv_compute_velocity(split1,dataset_short) 
scv_compute_velocity(split2,dataset_short) 

split1.write_h5ad(data_folder+'v4_test/v4_'+dataset_long+'/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_37-3_v4.h5ad')
split2.write_h5ad(data_folder+'v4_test/v4_'+dataset_long+'/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_37-7_v4.h5ad')
