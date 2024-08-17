import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import sctour as sct

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'sct'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 9815 × 2000
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split2_v2.h5ad') # 
#raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad") # 9815 × 53801

sc.pp.neighbors(split1, use_rep='X_TNODE', n_neighbors=15) # used 30 in first ery version
sc.tl.umap(split1)
sc.pp.neighbors(split2, use_rep='X_TNODE', n_neighbors=15) # used 30 in first ery version
sc.tl.umap(split2)
sct.vf.plot_vector_field(split1, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', show=False, size=100, alpha=0.2)
sct.vf.plot_vector_field(split2, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', show=False, size=100, alpha=0.2)

def sct_velocity_graph(adata):
    pos = adata.obsp['cosine_similarity'].copy()
    neg = adata.obsp['cosine_similarity'].copy()
    pos[pos < 0] = 0
    neg[neg > 0] = 0
    adata.uns['velocity_graph'] = pos
    adata.uns['velocity_graph_neg'] = neg

celltype_label = 'celltype'
sct_velocity_graph(split1)
sct_velocity_graph(split2)

def count_transition_matrix_pred_max(split1,split2,celltype_label,use_negative_cosines=False):
    Ncells = split1.shape[0]
    t1 = scv.utils.get_transition_matrix(split1,use_negative_cosines=use_negative_cosines)
    t2 = scv.utils.get_transition_matrix(split2,use_negative_cosines=use_negative_cosines)
    res_bool = []
    for cell_idx in range(Ncells):
        if cell_idx%1000==0: print(cell_idx)
        max_idx1 = np.where(t1[cell_idx].toarray()==np.max(t1[cell_idx]))[1][0]
        pred1 = split1.obs['celltype'].values[max_idx1]
        max_idx2 = np.where(t2[cell_idx].toarray()==np.max(t2[cell_idx]))[1][0]
        pred2 = split2.obs['celltype'].values[max_idx2]
        res_bool.append(pred1==pred2)
    return res_bool

res = count_transition_matrix_pred_max(split1,split2,celltype_label,use_negative_cosines=False)
np.sum(res), np.round(np.sum(res)/len(res),4) # (7029, 0.7161)

import collections
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))
"""
Blood progenitors 1: 0.5843
Blood progenitors 2: 0.678
Erythroid1: 0.5995
Erythroid2: 0.5542
Erythroid3: 0.9744
"""
for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))

