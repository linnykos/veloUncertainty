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

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

############## not revised yet
print_message_with_time('################## Read data') 
total = sc.read_h5ad(data_folder+dataset_short+'_total_allgenes.h5ad')
adata_split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')
gene_names = total.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

total.layers['test_spliced_original'] = total.layers['spliced'].copy()

S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()
S_mat_total = total.layers['spliced'].copy()
U_mat_total = total.layers['unspliced'].copy()

print_message_with_time('################## Run model on total')
scv_compute_velocity_ery(total) 
positions_total = [positions_dict[gene] for gene in total.var.index]
total.layers['spliced_original'] = S_mat_total[:,positions_total]
total.layers['unspliced_original'] = U_mat_total[:,positions_total]

print_message_with_time('################## Read model on split1')
adata_split1.layers['spliced_original'] = adata_split1.layers['spliced'].copy()
adata_split1.layers['unspliced_original'] = adata_split1.layers['unspliced'].copy()
scv_compute_velocity_ery(adata_split1)

print_message_with_time('################## Read model on split2')
adata_split2.layers['spliced_original'] = adata_split2.layers['spliced'].copy()
adata_split2.layers['unspliced_original'] = adata_split2.layers['unspliced'].copy()
scv_compute_velocity_ery(adata_split2)


raw = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad')
total.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
adata_split1.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
adata_split2.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()

# write data
print_message_with_time('################## Write data')
total.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_total_v4.h5ad')
adata_split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split1_v4.h5ad')
adata_split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split2_v4.h5ad')


scv.tl.velocity_confidence(total)
scv.tl.velocity_confidence(adata_split1)
scv.tl.velocity_confidence(adata_split2)
scv.tl.velocity_pseudotime(total)
scv.tl.velocity_pseudotime(adata_split1)
scv.tl.velocity_pseudotime(adata_split2)
scv.tl.latent_time(total)
scv.tl.latent_time(adata_split1)
scv.tl.latent_time(adata_split2)


total.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_total_v4_outputAdded.h5ad')
adata_split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split1_v4_outputAdded.h5ad')
adata_split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split2_v4_outputAdded.h5ad')

print_message_with_time('################## All done with the data')

### common genes being filtered out
common_genes = np.intersect1d(np.array(adata_split1.var.index), np.array(adata_split2.var.index)) 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0]))
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(adata_split1.var.index),np.array(total.var.index)).shape[0]))
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(adata_split2.var.index),np.array(total.var.index)).shape[0]))

plot_gene_correlation_between_splits(adata1=adata_split1,adata2=adata_split2,fig_folder=fig_folder,fig_path=dataset_short+'_'+method+'_corr_between_splits_colorSame.png')
