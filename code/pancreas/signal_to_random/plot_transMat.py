import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_transMat import *
import os
import re

dataset_short = 'pan'
dataset_long = 'pancreas'
celltype_label = 'celltype'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

## Mark
gene_set_name_prefix = 'Mark227'
gene_set_name = gene_set_name_prefix
split_seed = 317

for method_prefix in ['scv','utv','sct','velovi','velovi_woprep']:
    method = method_prefix + '_' + gene_set_name
    adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
    if 'scv' in method:
        path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total', s)][0]
    else: 
        path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total_*outputAdded', s)][0]
    total = sc.read_h5ad( path_total )
    celltypes = list(list(total.obs[celltype_label].cat.categories))
    trans_mat_mean, trans_mat_median = compute_celltype_transition_matrix_from_adata(adata_in=total, dataset=dataset_short)
    plot_transMat_heatmap(mat=trans_mat_mean, celltypes=celltypes, dataset=dataset_short, method=method, split_info='total',
                            prob_type='mean', fig_path=fig_folder+'seed'+str(split_seed)+"/"+method+"/",vmax=1)
    plot_transMat_heatmap(mat=trans_mat_median, celltypes=celltypes, dataset=dataset_short, method=method, split_info='total',
                            prob_type='median', fig_path=fig_folder+'seed'+str(split_seed)+"/"+method+"/",vmax=1)
    print(method+' done.')


## nMark
gene_set_name_prefix = 'nMark'
gene_set_name = gene_set_name_prefix

for method_prefix in ['scv','utv','sct','velovi','velovi_woprep']:
    for i in range(5):
        split_seed = [317,320,323,326,329][i]
        grid_seed = [227,230,233,236,239][i]
        method = method_prefix + '_' + gene_set_name + str(grid_seed)
        adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
        if 'scv' in method:
            path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total', s)][0]
        else: 
            path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total_*outputAdded', s)][0]
        total = sc.read_h5ad( path_total )
        celltypes = list(list(total.obs[celltype_label].cat.categories))
        trans_mat_mean, trans_mat_median = compute_celltype_transition_matrix_from_adata(adata_in=total, dataset=dataset_short)
        plot_transMat_heatmap(mat=trans_mat_mean, celltypes=celltypes, dataset=dataset_short, method=method, split_info='total',
                              prob_type='mean', fig_path=fig_folder+'seed'+str(split_seed)+"/"+method+"/",vmax=1)
        plot_transMat_heatmap(mat=trans_mat_median, celltypes=celltypes, dataset=dataset_short, method=method, split_info='total',
                              prob_type='median', fig_path=fig_folder+'seed'+str(split_seed)+"/"+method+"/",vmax=1)
        print(str(split_seed)+' done.')                                 
    print(method_prefix+' done.')
