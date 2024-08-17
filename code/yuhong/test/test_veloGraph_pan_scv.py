import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'scv'
dataset_long = 'pancreas'
dataset_short = 'pan'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 

"""
from scipy.stats import spearmanr
spearmanr(split1.uns['velocity_graph'][0].toarray(), split2.uns['velocity_graph'][0].toarray()).correlation

np.where(split1.uns['velocity_graph'][0].toarray()==np.max(split1.uns['velocity_graph'][0])) # 2890
np.where(split2.uns['velocity_graph'][0].toarray()==np.max(split2.uns['velocity_graph'][0])) # 1106

np.where(split1.uns['velocity_graph'][0].toarray()>0)[1]
"""

import collections

cell_idx=65
celltype_label = 'clusters'
res1 = []
res2 = []
#max1 = split1.obs[celltype_label][np.where(split1.uns['velocity_graph'][cell_idx].toarray()==np.max(split1.uns['velocity_graph'][cell_idx]))[1]][0]
#max2 = split2.obs[celltype_label][np.where(split2.uns['velocity_graph'][cell_idx].toarray()==np.max(split2.uns['velocity_graph'][cell_idx]))[1]][0]
for i in np.where(split1.uns['velocity_graph'][cell_idx].toarray()>0)[1]:
    celltype1=split1.obs[celltype_label][i]
    res1.append(celltype1)

for i in np.where(split2.uns['velocity_graph'][cell_idx].toarray()>0)[1]:
    celltype2=split2.obs[celltype_label][i]
    res2.append(celltype2)

counter1 = collections.Counter(res1)
counter2 = collections.Counter(res2)
counter1
counter2

counter1.most_common()[0][0]
counter2.most_common()[0][0]

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

count_velocity_graph_pred(split1,split2,celltype_label='clusters')
# pan+scv (3002, 664, 30)
# 3002/3696 = 0.8122
