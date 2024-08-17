#################################################################
## scv
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

def count_transition_matrix_pred(split1,split2,celltype_label):
    Ncells = split1.uns['velocity_graph'].shape[0]
    t1 = scv.utils.get_transition_matrix(split1)
    t2 = scv.utils.get_transition_matrix(split2)
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

res = count_transition_matrix_pred(split1,split2,celltype_label='celltype')
np.sum(res), np.round(np.sum(res)/len(res),4)
# ery+scv (7622, 0.7766)

celltype_label='celltype'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))
"""
Blood progenitors 1: 0.7255
Blood progenitors 2: 0.8667
Erythroid1: 0.7214
Erythroid2: 0.3852
Erythroid3: 0.9266
"""

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))


#################################################################
#################################################################
### utv
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import bbknn
import collections

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'utv'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 

def utv_compute_umap(adata):
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("************ batch correction done ************")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

utv_compute_umap(split1)
utv_compute_umap(split2)
utv_compute_umap(total)

def count_transition_matrix_pred(split1,split2,celltype_label):
    Ncells = split1.uns['velocity_graph'].shape[0]
    t1 = scv.utils.get_transition_matrix(split1)
    t2 = scv.utils.get_transition_matrix(split2)
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

res = count_transition_matrix_pred(split1,split2,celltype_label='celltype')
# ery+utv (8002, 1813, 0)=(Nsame, Ndiff, Nnan)
# 8002/9185 = 0.8153

celltype_label='celltype'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

np.array(res)[split1.obs[celltype_label].values==celltypes[0]]
np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[0]])/celltypes_counter[celltypes[0]]

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))
"""
Blood progenitors 1: 0.6613
Blood progenitors 2: 0.9033
Erythroid1: 0.7695
Erythroid2: 0.5072
Erythroid3: 0.9466
"""

#################################################################
#################################################################
### velovi
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

import matplotlib.pyplot as plt
import seaborn as sns
import collections

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_velovi import *

method = 'velovi'
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split1_v2.h5ad")
vae_split1 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_split1_v2.pt', split1)

split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split2_v2.h5ad")
vae_split2 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_split2_v2.pt', split2)

total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_total_v2.h5ad")
vae_total = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_total_v2.pt', total)

#raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")

## add velovi outputs to adata
add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)

## compute umap
compute_umap_ery(split1)
compute_umap_ery(split2)
compute_umap_ery(total)

import collections
def count_transition_matrix_pred(split1,split2,celltype_label):
    Ncells = split1.uns['velocity_graph'].shape[0]
    t1 = scv.utils.get_transition_matrix(split1)
    t2 = scv.utils.get_transition_matrix(split2)
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

res = count_transition_matrix_pred(split1,split2,celltype_label='celltype')
np.sum(res), np.round(np.sum(res)/len(res),4)
# ery+utv (6997, 0.7129)

celltype_label='celltype'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))

"""
Blood progenitors 1: 0.6212
Blood progenitors 2: 0.7171
Erythroid1: 0.7801
Erythroid2: 0.359
Erythroid3: 0.8024
"""

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))


#################################################################
#################################################################
### velovi_wopreprocess
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_velovi import *

method = 'velovi'
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/wopreprocess/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_v2.h5ad")
vae_split1 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_ery_velovi_split1_wopreprocess_v2.pt', split1)

split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_v2.h5ad")
vae_split2 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_ery_velovi_split2_wopreprocess_v2.pt', split2)

total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_v2.h5ad")
vae_total = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_ery_velovi_total_wopreprocess_v2.pt', total)

#raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")

add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)

compute_umap_ery(split1)
compute_umap_ery(split2)
compute_umap_ery(total)

split1.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_outputAdded_v2.h5ad")
split2.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_outputAdded_v2.h5ad")
total.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_outputAdded_v2.h5ad")

import collections
def count_transition_matrix_pred(split1,split2,celltype_label):
    Ncells = split1.uns['velocity_graph'].shape[0]
    t1 = scv.utils.get_transition_matrix(split1)
    t2 = scv.utils.get_transition_matrix(split2)
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

res = count_transition_matrix_pred(split1,split2,celltype_label='celltype')
np.sum(res), np.round(np.sum(res)/len(res),4)
# ery+utv (7623, 0.7767)

celltype_label='celltype'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))

"""
Blood progenitors 1: 0.7255
Blood progenitors 2: 0.8667
Erythroid1: 0.7214
Erythroid2: 0.3861
Erythroid3: 0.9266
"""

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))




