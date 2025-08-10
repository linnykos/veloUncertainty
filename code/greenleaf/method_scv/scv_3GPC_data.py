import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'scv_GPC'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

## different from pancreas read adata function
print('################## Read data') 
total = sc.read(data_folder+'glf_total_GPC.h5ad')
split1 = sc.read(data_folder+'seed317_glf_split1_GPC.h5ad')
split2 = sc.read(data_folder+'seed317_glf_split2_GPC.h5ad')

split1.layers['spliced_original'] = split1.layers['spliced'].copy()
split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()
split2.layers['spliced_original'] = split1.layers['spliced'].copy()
split2.layers['unspliced_original'] = split1.layers['unspliced'].copy()
total.layers['spliced_original'] = total.layers['spliced'].copy()
total.layers['unspliced_original'] = total.layers['unspliced'].copy()

print('################## start velocity') 
scv_compute_velocity(total,dataset_short) 
scv_compute_velocity(split1,dataset_short) 
scv_compute_velocity(split2,dataset_short) 
print('################## end velocity') 

total.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_total_GPC.h5ad')
split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')

