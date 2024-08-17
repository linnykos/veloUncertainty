import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import collections

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

dataset_short = 'ery'
dataset_long = 'erythroid'
method = 'scv'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/scv/" 

total = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_total_v2.h5ad')
split1 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split1_v2.h5ad')
split2 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split2_v2.h5ad')

def count_velocity_graph_pred(split1,split2,celltype_label):
    Ncells = split1.uns['velocity_graph'].shape[0]
    Nsame = 0
    Ndiff = 0
    Nnan = 0
    for cell_idx in range(Ncells):
        print(cell_idx)
        res1 = []
        for i in np.where(split1.uns['velocity_graph'][cell_idx].toarray()>0)[1]:
            celltype1=split1.obs[celltype_label][i]
            res1.append(celltype1)
        if len(res1)==0: 
            Nnan += 1
            continue
        res2 = []
        for i in np.where(split2.uns['velocity_graph'][cell_idx].toarray()>0)[1]:
            celltype2=split2.obs[celltype_label][i]
            res2.append(celltype2)
        if len(res2)==0: 
            Nnan += 1
            continue
        counter1 = collections.Counter(res1)
        counter2 = collections.Counter(res2)
        pred1 = counter1.most_common()[0][0]
        pred2 = counter2.most_common()[0][0]
        if pred1==pred2: Nsame += 1
        else: Ndiff += 1
    return Nsame, Ndiff, Nnan

count_velocity_graph_pred(split1,split2,celltype_label='celltype')
# ery+scv (7313, 2345, 157)
# 7313/9815 = 0.7451

t1 = scv.utils.get_transition_matrix(split1)
np.where(t1[0].toarray()>0)[1]

cell_idx = 103
np.where(t1[cell_idx].toarray()==np.max(t1[cell_idx].todense()))[1]
np.where(split1.uns['velocity_graph'][cell_idx].todense()==np.max(split1.uns['velocity_graph'][cell_idx].toarray()))[1]


np.where(split1.uns['velocity_graph'][0].toarray()>0)[1]

def count_transition_matrix_pred(split1,split2,celltype_label):
    Ncells = split1.uns['velocity_graph'].shape[0]
    t1 = scv.utils.get_transition_matrix(split1)
    t2 = scv.utils.get_transition_matrix(split2)
    Nsame = 0
    Ndiff = 0
    Nnan = 0
    for cell_idx in range(Ncells):
        print(cell_idx)
        res1 = []
        for i in np.where(t1[cell_idx].toarray()>0)[1]:
            celltype1=split1.obs[celltype_label][i]
            res1.append(celltype1)
        if len(res1)==0: 
            Nnan += 1
            continue
        res2 = []
        for i in np.where(t2[cell_idx].toarray()>0)[1]:
            celltype2=split2.obs[celltype_label][i]
            res2.append(celltype2)
        if len(res2)==0: 
            Nnan += 1
            continue
        counter1 = collections.Counter(res1)
        counter2 = collections.Counter(res2)
        pred1 = counter1.most_common()[0][0]
        pred2 = counter2.most_common()[0][0]
        if pred1==pred2: Nsame += 1
        else: Ndiff += 1
    return Nsame, Ndiff, Nnan

count_transition_matrix_pred(split1,split2,celltype_label='celltype')
# ery+scv (7622, 2193, 0)
# 0.7766
