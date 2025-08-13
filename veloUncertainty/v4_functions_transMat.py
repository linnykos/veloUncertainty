import scvelo as scv
import numpy as np
import collections

import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
#from v4_functions import get_umap_sct

## transMat
import seaborn as sns
def compute_celltype_idx_dict(adata, celltype_label):
    celltype_idx = {}
    celltypes = adata.obs[celltype_label].cat.categories
    for celltype_i in range(len(celltypes)):
        celltype = celltypes[celltype_i]
        cell_idx = np.where([adata.obs[celltype_label][i]==celltype for i in range(adata.shape[0])])[0]
        celltype_idx[celltype] = list(cell_idx)
    return celltype_idx

def compute_celltype_transition_matrix_from_vk(vk, celltype_label, celltype_idx):
    celltypes = list(celltype_idx.keys())
    trans_prob_median_mat = np.zeros((len(celltypes), len(celltypes)))
    trans_prob_mean_mat = np.zeros((len(celltypes), len(celltypes)))
    for celltype_idx_i in range(len(celltypes)):
        for celltype_idx_j in range(len(celltypes)):
            celltype_i = celltypes[celltype_idx_i]
            celltype_j = celltypes[celltype_idx_j]
            celltype_mat_ij = vk.transition_matrix[celltype_idx[celltype_i],][:,celltype_idx[celltype_j]]
            probs = [np.sum(celltype_mat_ij[i,]) for i in range(celltype_mat_ij.shape[0])]
            trans_prob_median_mat[celltype_idx_i,celltype_idx_j] = np.round(np.median(probs),4)
            trans_prob_mean_mat[celltype_idx_i,celltype_idx_j] = np.round(np.mean(probs),4)
    return trans_prob_mean_mat, trans_prob_median_mat

def compute_celltype_transition_matrix_from_adata(adata_in, dataset, copy=True, celltype_label=None):
    adata = None
    if copy: adata = adata_in
    else: adata = adata_in
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()
    if celltype_label==None:
        if 'pan' in dataset: celltype_label = 'clusters'
        elif 'ery' in dataset: celltype_label = 'celltype'
        elif 'larry' in dataset: celltype_label = 'state_info'
    celltype_idx = compute_celltype_idx_dict(adata, celltype_label)
    return compute_celltype_transition_matrix_from_vk(vk, celltype_label, celltype_idx)

def plot_transMat_heatmap(mat, celltypes, dataset, method, split_info, fig_path, prob_type,cmap='Blues',vmin=0,vmax=None):
    plt.figure(figsize=(10, 8), facecolor='w')
    sns.heatmap(mat,annot=True,fmt=".2f",cmap=cmap,cbar=True,square=True,vmax=vmax,vmin=vmin,xticklabels=celltypes,yticklabels=celltypes)
    plt.tick_params(colors='black')
    plt.xticks(rotation=45, ha='right')
    plt.title('transition matrix, prob '+prob_type+', '+dataset+'+'+method+' '+split_info)
    plt.savefig(fig_path+'transMat/'+dataset+'_'+method+'_'+prob_type+'_'+split_info+'.png')
    plt.clf()

def plot_transMat_heatmap_from_adata(split1, split2, total, method, dataset_short, fig_folder,celltype_label=None):
    if ('pan' in dataset_short): celltype_label = 'clusters'
    elif ('ery' in dataset_short): celltype_label = 'celltype'
    elif ('larry' in dataset_short): celltype_label = 'state_info'
    celltypes = list(list(split1.obs[celltype_label].cat.categories))
    if (method=='sct'):
        get_umap_sct(split1)
        get_umap_sct(split2)
        get_umap_sct(total)
    #celltype_idx = compute_celltype_idx_dict(total, celltype_label)
    trans_mat1_mean, trans_mat1_median = compute_celltype_transition_matrix_from_adata(adata_in=split1, dataset=dataset_short)
    trans_mat2_mean, trans_mat2_median = compute_celltype_transition_matrix_from_adata(adata_in=split2, dataset=dataset_short)
    trans_mat_mean, trans_mat_median = compute_celltype_transition_matrix_from_adata(adata_in=total, dataset=dataset_short)
    # mean
    plot_transMat_heatmap(mat=trans_mat1_mean, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='split1',prob_type='mean', fig_path=fig_folder,vmax=1)
    plot_transMat_heatmap(mat=trans_mat2_mean, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='split2',prob_type='mean', fig_path=fig_folder,vmax=1)
    plot_transMat_heatmap(mat=np.abs(trans_mat1_mean-trans_mat2_mean), celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='diff_abs',prob_type='mean', fig_path=fig_folder, cmap='Reds',vmax=.3)
    plot_transMat_heatmap(mat=trans_mat_mean, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='total',prob_type='mean', fig_path=fig_folder,vmax=1)
    # median
    plot_transMat_heatmap(mat=trans_mat1_median, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='split1',prob_type='median', fig_path=fig_folder,vmax=1)
    plot_transMat_heatmap(mat=trans_mat2_median, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='split2',prob_type='median', fig_path=fig_folder,vmax=1)
    plot_transMat_heatmap(mat=np.abs(trans_mat1_median-trans_mat2_median), celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='diff_abs',prob_type='median', fig_path=fig_folder, cmap='Reds',vmax=.3)
    plot_transMat_heatmap(mat=trans_mat_median, celltypes=celltypes, dataset=dataset_short, 
                     method=method, split_info='total',prob_type='median', fig_path=fig_folder,vmax=1)
