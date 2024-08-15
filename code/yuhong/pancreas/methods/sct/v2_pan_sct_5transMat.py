import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_transMat import *

# the all-in-one version
method = 'sct'
dataset_long = 'pancreas'
dataset_short = 'pan'
celltype_label = 'clusters'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
plot_transMat_heatmap_from_adata(split1=split1, split2=split2, total=total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)


exit()
##########################################
# specific steps
celltypes = list(list(split1.obs[celltype_label].cat.categories))

get_umap_sct(split1)
get_umap_sct(split2)
get_umap_sct(total)

celltype_idx = compute_celltype_idx_dict(total, celltype_label)
trans_mat1_mean, trans_mat1_median = compute_celltype_transition_matrix_from_adata(adata_in=split1, dataset=dataset_short)
trans_mat2_mean, trans_mat2_median = compute_celltype_transition_matrix_from_adata(adata_in=split2, dataset=dataset_short)
trans_mat_mean, trans_mat_median = compute_celltype_transition_matrix_from_adata(adata_in=total, dataset=dataset_short)


plot_transMat_heatmap(mat=trans_mat1_mean, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='split1',prob_type='mean', fig_path=fig_folder)
plot_transMat_heatmap(mat=trans_mat2_mean, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='split2',prob_type='mean', fig_path=fig_folder)
plot_transMat_heatmap(mat=np.abs(trans_mat1_mean-trans_mat2_mean), celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='diff_abs',prob_type='mean', fig_path=fig_folder)

plot_transMat_heatmap(mat=trans_mat_mean, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='total',prob_type='mean', fig_path=fig_folder)


