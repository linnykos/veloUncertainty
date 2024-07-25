import scvelo as scv
import numpy as np
import collections

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
