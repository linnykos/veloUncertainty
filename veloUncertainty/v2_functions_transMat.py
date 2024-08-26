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
from v2_functions import get_umap_sct

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


"""
def print_alignment_by_celltype(adata,celltype_label,res_pred,c_same_neighbor,print_celltype=True):
    celltypes = adata.obs[celltype_label].cat.categories
    celltypes_counter = collections.Counter(adata.obs[celltype_label])
    for ct in celltypes:
        idx = np.where(celltypes==ct)[0][0]
        ct_idx_set = adata.obs[celltype_label].values==celltypes[idx]
        num = np.sum(np.array(res_pred)[ct_idx_set]) - np.sum(np.array(c_same_neighbor)[ct_idx_set])
        denom = celltypes_counter[celltypes[idx]] - np.sum(np.array(c_same_neighbor)[ct_idx_set])
        if print_celltype: print(ct+': '+str(np.round(num / denom, 4)))
        else: print(str(np.round(num / denom, 4)))

def calculate_transMat_alignment(split1,split2,celltype_label,use_negative_cosines=False,print_celltype=True,return_array=False,correct_c_same_neighbor=True):
    Ncells = split1.shape[0]
    t1 = scv.utils.get_transition_matrix(split1,use_negative_cosines=use_negative_cosines)
    t2 = scv.utils.get_transition_matrix(split2,use_negative_cosines=use_negative_cosines)
    res_pred = []
    res_c = []
    for cell_idx in range(Ncells):
        if cell_idx%1000==0: print(cell_idx)
        # count_same_neighbors
        nbr1_idx = np.where(split1.obsp['connectivities'][cell_idx].todense()!=0)[1]
        nbr2_idx = np.where(split2.obsp['connectivities'][cell_idx].todense()!=0)[1]
        nbr1_counter = collections.Counter(split1.obs[celltype_label][nbr1_idx])
        nbr2_counter = collections.Counter(split2.obs[celltype_label][nbr2_idx])
        cur_c = ( (len(nbr1_counter)==1 and len(nbr2_counter)==1 and (list(nbr1_counter)[0] == list(nbr1_counter)[0])) )
        res_c.append(cur_c)
        # count_same_pred_by_transition_matrix_max
        max_idx1 = np.where(t1[cell_idx].toarray()==np.max(t1[cell_idx]))[1][0]
        pred1 = split1.obs[celltype_label].values[max_idx1]
        max_idx2 = np.where(t2[cell_idx].toarray()==np.max(t2[cell_idx]))[1][0]
        pred2 = split2.obs[celltype_label].values[max_idx2]
        res_pred.append(pred1==pred2)
    if correct_c_same_neighbor==False: res_c = [0]*len(res_c)
    c = np.sum(res_c)
    print('### Results:')
    print('Number of same celltype prediction = '+str(np.sum(res_pred)-c))
    print('Proportion of same celltype prediction = '+str(np.round( (np.sum(res_pred)-c)/(len(res_pred)-c), 4 ) ))
    print_alignment_by_celltype(adata=split1,celltype_label=celltype_label,res_pred=res_pred,c_same_neighbor=res_c,print_celltype=print_celltype)
    if return_array: return res_c, res_pred
    else: return
"""