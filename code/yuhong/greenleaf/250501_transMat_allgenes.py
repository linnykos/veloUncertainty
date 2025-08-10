import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_transMat import *
import os
import re

dataset_short = 'glf'
dataset_long = 'greenleaf'
celltype_label = 'cluster_name'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

## allgenes
split_seed = 317

for method in ['scv','utv','sct','velovi','velovi_woprep']:
for method in ['utv','sct','velovi','velovi_woprep']:
    adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
    if ('scv' in method) or ('utv' in method):
        path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total', s)][0]
    else: 
        path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total_.*outputAdded', s)][0]
    total = sc.read_h5ad( path_total )
    celltypes = list(list(total.obs[celltype_label].cat.categories))
    trans_mat_mean, trans_mat_median = compute_celltype_transition_matrix_from_adata(adata_in=total, dataset=dataset_short, celltype_label=celltype_label)
    plot_transMat_heatmap(mat=trans_mat_mean, celltypes=celltypes, dataset=dataset_short, method=method, split_info='total',
                            prob_type='mean', fig_path=fig_folder+'seed'+str(split_seed)+"/"+method+"/",vmax=1)
    print(method+' done.')



def compute_Ngenes_intersect(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    if "scv" in method:
        velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
        velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0])) 
    #return common_genes_velocity.shape[0]

def compute_Ngenes_union(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    velo_split1 = pd.DataFrame(adata_split1.layers['velocity'], columns=velo_genes_split1)
    velo_split2 = pd.DataFrame(adata_split2.layers['velocity'], columns=velo_genes_split2)
    if 'scv' in method:
        velo_genes_split1 = velo_genes_split1[~np.isnan(velo_split1.loc[0])] #adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
        velo_genes_split2 = velo_genes_split2[~np.isnan(velo_split2.loc[0])] #adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    union_genes_velo = np.union1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Size of the union of genes for velocity computation in splits = '+str(union_genes_velo.shape[0])) 
    
    #return union_genes_velo.shape[0]


## GPC
gene_set_name_prefix = 'GPC'
gene_set_name = gene_set_name_prefix
split_seed = 317

for method_prefix in ['scv','utv','sct','velovi','velovi_woprep']:
    method = method_prefix + '_' + gene_set_name
    print(method)
    adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
    if ('scv' in method) or ('utv' in method):
        path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1', s)][0]
        path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2', s)][0]
        #path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total', s)][0]
    else: 
        path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1_.*outputAdded.*h5ad', s)][0]
        path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2_.*outputAdded.*h5ad', s)][0]
        #path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total_.*outputAdded', s)][0]
    #total = sc.read_h5ad( path_total )
    split1 = sc.read_h5ad( path_split1 )
    split2 = sc.read_h5ad( path_split2 )
    celltypes = list(list(total.obs[celltype_label].cat.categories))
    compute_Ngenes_intersect(split1,split2,method)
    compute_Ngenes_union(split1,split2,method)
    print(method+' done.')

"""
GPC
scv_GPC
Number of overlapped genes for velocity computation in splits = 67
Size of the union of genes for velocity computation in splits = 77
scv_GPC done.
utv_GPC
Number of overlapped genes for velocity computation in splits = 144
Size of the union of genes for velocity computation in splits = 150
utv_GPC done.
sct_GPC
Number of overlapped genes for velocity computation in splits = 184
Size of the union of genes for velocity computation in splits = 184
sct_GPC done.
velovi_GPC
Number of overlapped genes for velocity computation in splits = 65
Size of the union of genes for velocity computation in splits = 78
velovi_GPC done.
velovi_woprep_GPC
Number of overlapped genes for velocity computation in splits = 144
Size of the union of genes for velocity computation in splits = 150
velovi_woprep_GPC done.
"""

## nGPC
gene_set_name_prefix = 'nGPCblk'
grid_seed = 227
gene_set_name = gene_set_name_prefix + str(grid_seed)
split_seed = 317

gene_set_name_prefix = 'nGPCblk'
grid_seed = 233
gene_set_name = gene_set_name_prefix + str(grid_seed)
split_seed = 323

for method_prefix in ['scv','utv','sct','velovi','velovi_woprep']:
    method = method_prefix + '_' + gene_set_name
    print(method)
    adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
    if ('scv' in method) or ('utv' in method):
        path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1', s)][0]
        path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2', s)][0]
        #path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total', s)][0]
    else: 
        path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1_.*outputAdded.*h5ad', s)][0]
        path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2_.*outputAdded.*h5ad', s)][0]
        #path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total_.*outputAdded', s)][0]
    #total = sc.read_h5ad( path_total )
    split1 = sc.read_h5ad( path_split1 )
    split2 = sc.read_h5ad( path_split2 )
    celltypes = list(list(total.obs[celltype_label].cat.categories))
    compute_Ngenes_intersect(split1,split2,method)
    compute_Ngenes_union(split1,split2,method)
    print(method+' done.')

"""
scv_nGPCblk227
Number of overlapped genes for velocity computation in splits = 72
Size of the union of genes for velocity computation in splits = 75
scv_nGPCblk227 done.
utv_nGPCblk227
Number of overlapped genes for velocity computation in splits = 219
Size of the union of genes for velocity computation in splits = 224
utv_nGPCblk227 done.
sct_nGPCblk227
Number of overlapped genes for velocity computation in splits = 368
Size of the union of genes for velocity computation in splits = 368
sct_nGPCblk227 done.
velovi_nGPCblk227
Number of overlapped genes for velocity computation in splits = 68
Size of the union of genes for velocity computation in splits = 72
velovi_nGPCblk227 done.
velovi_woprep_nGPCblk227
Number of overlapped genes for velocity computation in splits = 219
Size of the union of genes for velocity computation in splits = 224
velovi_woprep_nGPCblk227 done.
"""

"""
scv_nGPCblk230
Number of overlapped genes for velocity computation in splits = 58
Size of the union of genes for velocity computation in splits = 66
scv_nGPCblk230 done.
utv_nGPCblk230
Number of overlapped genes for velocity computation in splits = 219
Size of the union of genes for velocity computation in splits = 226
utv_nGPCblk230 done.
sct_nGPCblk230
Number of overlapped genes for velocity computation in splits = 368
Size of the union of genes for velocity computation in splits = 368
sct_nGPCblk230 done.
velovi_nGPCblk230
Number of overlapped genes for velocity computation in splits = 56
Size of the union of genes for velocity computation in splits = 66
velovi_nGPCblk230 done.
velovi_woprep_nGPCblk230
Number of overlapped genes for velocity computation in splits = 219
Size of the union of genes for velocity computation in splits = 226
velovi_woprep_nGPCblk230 done.
"""

"""
scv_nGPCblk233
Number of overlapped genes for velocity computation in splits = 53
Size of the union of genes for velocity computation in splits = 62
scv_nGPCblk233 done.
utv_nGPCblk233
Number of overlapped genes for velocity computation in splits = 215
Size of the union of genes for velocity computation in splits = 220
utv_nGPCblk233 done.
sct_nGPCblk233
Number of overlapped genes for velocity computation in splits = 368
Size of the union of genes for velocity computation in splits = 368
sct_nGPCblk233 done.
velovi_nGPCblk233
Number of overlapped genes for velocity computation in splits = 51
Size of the union of genes for velocity computation in splits = 60
velovi_nGPCblk233 done.
velovi_woprep_nGPCblk233
Number of overlapped genes for velocity computation in splits = 215
Size of the union of genes for velocity computation in splits = 220
velovi_woprep_nGPCblk233 done.
"""
