import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'scv'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

## different from pancreas read adata function
print_message_with_time('################## Read data') 
total = sc.read_h5ad(data_folder+dataset_short+'_total_allgenes.h5ad')
split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')


split1.layers['spliced_original'] = split1.layers['spliced'].copy()
split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()
split2.layers['spliced_original'] = split1.layers['spliced'].copy()
split2.layers['unspliced_original'] = split1.layers['unspliced'].copy()
total.layers['spliced_original'] = total.layers['spliced'].copy()
total.layers['unspliced_original'] = total.layers['unspliced'].copy()

print_message_with_time('################## start velocity') 
scv_compute_velocity(total,dataset_short) 
scv_compute_velocity(split1,dataset_short) 
scv_compute_velocity(split2,dataset_short) 
print_message_with_time('################## end velocity') 


total.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_total_v4.h5ad')
split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')

######################


