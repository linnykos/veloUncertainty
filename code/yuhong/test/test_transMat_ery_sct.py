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
raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad") # 9815 × 53801

"""
AnnData object with n_obs x n_vars = 9815 x 2000
    obs: 'sample', 'stage', 'sequencing.batch', 'theiler', 'celltype', 'n_genes_by_counts', 'total_counts', 'ptime'
    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand', 'MURK_gene', 'Δm', 'scaled Δm', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'
    uns: 'celltype_colors', 'hvg'
    obsm: 'X_TNODE', 'X_VF', 'X_pcaOriginal', 'X_umapOriginal'
    layers: 'spliced', 'unspliced', 'velocity'
"""
t1 = sct.vf.cosine_similarity(split1,'X_TNODE')
"""
AnnData object with n_obs x n_vars = 9815 x 2000
    obs: 'sample', 'stage', 'sequencing.batch', 'theiler', 'celltype', 'n_genes_by_counts', 'total_counts', 'ptime'
    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand', 'MURK_gene', 'Δm', 'scaled Δm', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'
    uns: 'celltype_colors', 'hvg', 'neighbors'
    obsm: 'X_TNODE', 'X_VF', 'X_pcaOriginal', 'X_umapOriginal'
    layers: 'spliced', 'unspliced', 'velocity'
    obsp: 'distances', 'connectivities'
"""

t2 = sct.vf.cosine_similarity(split2,'X_TNODE')

sct.vf.plot_vector_field(split1, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', show=False, size=100, alpha=0.2)
sct.vf.plot_vector_field(split2, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', show=False, size=100, alpha=0.2)
"""
AnnData object with n_obs × n_vars = 9815 × 2000
    obs: 'sample', 'stage', 'sequencing.batch', 'theiler', 'celltype', 'n_genes_by_counts', 'total_counts', 'ptime'
    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand', 'MURK_gene', 'Δm', 'scaled Δm', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'
    uns: 'celltype_colors', 'hvg', 'neighbors'
    obsm: 'X_TNODE', 'X_VF', 'X_pcaOriginal', 'X_umapOriginal', 'X_umap', 'X_DV'
    layers: 'spliced', 'unspliced', 'velocity'
    obsp: 'distances', 'connectivities', 'cosine_similarity'
"""
split1 = split1[np.argsort(split1.obs['ptime'].values), :]
sc.pp.neighbors(split1, use_rep='X_TNODE', n_neighbors=15) # used 30 in first ery version
sc.tl.umap(split1)

split1.obsp['cosine_similarity']
split2.obsp['cosine_similarity']

np.sum(split1.obsp['cosine_similarity'][0]>0)
np.sum(split2.obsp['cosine_similarity'][0]>0)

# https://github.com/LiQian-XC/sctour/blob/86e21fe356012776d55a07f7c2fa7013dc19def8/sctour/vector_field.py#L15

def count_sct_cosine_similarity_pred_pos(split1,split2,celltype_label):
    import collections
    Ncells = split1.shape[0]
    t1 = split1.obsp['cosine_similarity']
    t2 = split2.obsp['cosine_similarity']
    res_bool = []
    for cell_idx in range(Ncells):
        if cell_idx%1000==0: print(cell_idx)
        res1 = []
        for i in np.where(t1[cell_idx].toarray()>0)[1]:
            celltype1=split1.obs[celltype_label][i]
            res1.append(celltype1)
            counter1 = collections.Counter(res1)
        res2 = []
        for i in np.where(t2[cell_idx].toarray()>0)[1]:
            celltype2=split2.obs[celltype_label][i]
            res2.append(celltype2)
            counter2 = collections.Counter(res2)
        pred1 = counter1.most_common()[0][0]
        pred2 = counter2.most_common()[0][0]
        res_bool.append(pred1==pred2)
    return res_bool

res = count_sct_cosine_similarity_pred_pos(split1,split2,celltype_label='celltype')
np.sum(res), np.round(np.sum(res)/len(res),4)
# (7760, 0.7906)

def count_sct_cosine_similarity_pred_neg(split1,split2,celltype_label,use_negative_cosines=False):
    import collections
    Ncells = split1.shape[0]
    t1 = split1.obsp['cosine_similarity']
    t2 = split2.obsp['cosine_similarity']
    res_bool = []
    for cell_idx in range(Ncells):
        if cell_idx%1000==0: print(cell_idx)
        res1 = []
        for i in np.where(t1[cell_idx].toarray()<0)[1]:
            celltype1=split1.obs[celltype_label][i]
            res1.append(celltype1)
            counter1 = collections.Counter(res1)
        res2 = []
        for i in np.where(t2[cell_idx].toarray()<0)[1]:
            celltype2=split2.obs[celltype_label][i]
            res2.append(celltype2)
            counter2 = collections.Counter(res2)
        pred1 = counter1.most_common()[0][0]
        pred2 = counter2.most_common()[0][0]
        res_bool.append(pred1==pred2)
    return res_bool

res_neg = count_sct_cosine_similarity_pred_neg(split1,split2,celltype_label='celltype')
np.sum(res_neg), np.round(np.sum(res_neg)/len(res_neg),4)
# (7553, 0.7695)

split1.obsp['cosine_similarity'][split1.obsp['cosine_similarity']>0]

pos1 = split1.obsp['cosine_similarity'].copy()
neg1 = split1.obsp['cosine_similarity'].copy()
pos1[pos1 < 0] = 0
neg1[neg1 > 0] = 0
split1.uns['velocity_graph'] = pos1
split1.uns['velocity_graph_neg'] = neg1
t1 = scv.utils.get_transition_matrix(split1)

pos2 = split2.obsp['cosine_similarity'].copy()
neg2 = split2.obsp['cosine_similarity'].copy()
pos2[pos2 < 0] = 0
neg2[neg2 > 0] = 0
split2.uns['velocity_graph'] = pos2
split2.uns['velocity_graph_neg'] = neg2
t2 = scv.utils.get_transition_matrix(split2)

import collections
def count_transition_matrix_pred(split1,split2,celltype_label,use_negative_cosines=False):
    Ncells = split1.uns['velocity_graph'].shape[0]
    t1 = scv.utils.get_transition_matrix(split1,use_negative_cosines=use_negative_cosines)
    t2 = scv.utils.get_transition_matrix(split2,use_negative_cosines=use_negative_cosines)
    res_bool = []
    for cell_idx in range(Ncells):
        if cell_idx%1000==0: print(cell_idx)
        res1 = []
        for i in np.where(t1[cell_idx].toarray()>0)[1]:
            celltype1=split1.obs[celltype_label][i]
            res1.append(celltype1)
            counter1 = collections.Counter(res1)
        res2 = []
        for i in np.where(t2[cell_idx].toarray()>0)[1]:
            celltype2=split2.obs[celltype_label][i]
            res2.append(celltype2)
            counter2 = collections.Counter(res2)
        pred1 = counter1.most_common()[0][0]
        pred2 = counter2.most_common()[0][0]
        res_bool.append(pred1==pred2)
    return res_bool


import collections
res = count_transition_matrix_pred(split1,split2,celltype_label='celltype',use_negative_cosines=True)
np.sum(res), np.round(np.sum(res)/len(res),4)
# use_negative_cosines=False: (7801, 0.7948)
# use_negative_cosines=True: (7790, 0.7937)

celltype_label='celltype'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))

"""
use_negative_cosines=False
Blood progenitors 1: 0.7737
Blood progenitors 2: 0.7683
Erythroid1: 0.8255
Erythroid2: 0.4638
Erythroid3: 0.9262

use_negative_cosines=True
Blood progenitors 1: 0.7079
Blood progenitors 2: 0.7528
Erythroid1: 0.7812
Erythroid2: 0.547
Erythroid3: 0.9655
"""