import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'scv'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_test/v4_'+dataset_long+'/'

############## not revised yet
print_message_with_time('################## Read data') 
adata_split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_allgenes_37-3.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_allgenes_37-7.h5ad')
gene_names = adata_split1.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()

print_message_with_time('################## Read model on split1')
adata_split1.layers['spliced_original'] = adata_split1.layers['spliced'].copy()
adata_split1.layers['unspliced_original'] = adata_split1.layers['unspliced'].copy()
scv_compute_velocity_ery(adata_split1)

print_message_with_time('################## Read model on split2')
adata_split2.layers['spliced_original'] = adata_split2.layers['spliced'].copy()
adata_split2.layers['unspliced_original'] = adata_split2.layers['unspliced'].copy()
scv_compute_velocity_ery(adata_split2)


raw = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad')
adata_split1.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
adata_split2.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()

# write data
print_message_with_time('################## Write data')
adata_split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_37-3_v4.h5ad')
adata_split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_37-7_v4.h5ad')

print_message_with_time('################## All done with the data')


