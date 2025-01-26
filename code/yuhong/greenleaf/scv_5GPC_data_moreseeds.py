import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'scv_GPC'

## different from pancreas read adata function
for split_seed in [320,323,326,329]:
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    print('################## Read data') 
    split1 = sc.read(data_folder+'seed'+str(split_seed)+'_glf_split1_GPC.h5ad')
    split2 = sc.read(data_folder+'seed'+str(split_seed)+'_glf_split2_GPC.h5ad')
    split1.layers['spliced_original'] = split1.layers['spliced'].copy()
    split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()
    split2.layers['spliced_original'] = split1.layers['spliced'].copy()
    split2.layers['unspliced_original'] = split1.layers['unspliced'].copy()
    print('################## start velocity') 
    scv_compute_velocity(split1,dataset_short) 
    scv_compute_velocity(split2,dataset_short) 
    print('################## end velocity') 
    split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
    split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')

