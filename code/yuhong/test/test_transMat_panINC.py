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

method = 'scv'
dataset_long = 'pancreasINC'
dataset_short = 'panINC'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
raw = scv.read(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")

cell_index = np.array(np.where(raw.obs['clusters']!='Pre-endocrine')[0])
def create_adata_INC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}#adata_old.uns['clusters_colors'].copy()
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

S_raw = raw.layers['spliced'][cell_index,:]
U_raw = raw.layers['unspliced'][cell_index,:]
raw = create_adata_INC(S=S_raw,U=U_raw,adata_old=raw)

total.uns['clusters_colors'] = split1.uns['clusters_colors'].copy()

##################################
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

res = count_transition_matrix_pred(split1,split2,celltype_label='clusters')
# panINC+scv, Nsame=2808
# 2808/3104 = 0.9046

celltype_label='clusters'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))
"""
Ductal: 0.9978
Ngn3 low EP: 0.9618
Ngn3 high EP: 0.9486
Beta: 0.9679
Alpha: 0.7422
Delta: 0.4286
Epsilon: 0.5211
"""


#################################################
#### utv
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

method = 'utv'
dataset_long = 'pancreasINC'
dataset_short = 'panINC'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
raw = scv.read(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")

cell_index = np.array(np.where(raw.obs['clusters']!='Pre-endocrine')[0])
def create_adata_INC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

S_raw = raw.layers['spliced'][cell_index,:]
U_raw = raw.layers['unspliced'][cell_index,:]
raw = create_adata_INC(S=S_raw,U=U_raw,adata_old=raw)

total.uns['clusters_colors'] = split1.uns['clusters_colors'].copy()

## compute umap
sc.tl.umap(total)
sc.tl.umap(split1)
sc.tl.umap(split2)

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

res = count_transition_matrix_pred(split1,split2,celltype_label='clusters')
np.sum(res), np.sum(res)/len(res)
# panINC+utv, Nsame=2841
# 2841/3104 = 0.9153

celltype_label='clusters'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))
"""
Ductal: 1.0
Ngn3 low EP: 1.0
Ngn3 high EP: 0.9564
Beta: 0.9695
Alpha: 0.7193
Delta: 0.8571
Epsilon: 0.493
"""

#################################################
#### velovi
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import collections

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_velovi import *

method = 'velovi'
dataset_short = 'panINC'
dataset_long = 'pancreasINC'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split1_v2.h5ad")
vae_split1 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_panINC_velovi_split1_v2.pt', split1)

split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split2_v2.h5ad")
vae_split2 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_panINC_velovi_split2_v2.pt', split2)

total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_total_v2.h5ad")
vae_total = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_panINC_velovi_total_v2.pt', total)

raw = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
cell_index = np.array(np.where(raw.obs['clusters']!='Pre-endocrine')[0])
def create_adata_INC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

S_raw = raw.layers['spliced'][cell_index,:]
U_raw = raw.layers['unspliced'][cell_index,:]
raw = create_adata_INC(S=S_raw,U=U_raw,adata_old=raw)

total.uns['clusters_colors'] = split1.uns['clusters_colors'].copy()

## add velovi outputs to adata
add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)

## compute umap
compute_umap_pan(split1)
compute_umap_pan(split2)
compute_umap_pan(total)

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

res = count_transition_matrix_pred(split1,split2,celltype_label='clusters')
np.sum(res), np.sum(res)/len(res)
# panINC+velovi, Nsame=2670
# /3104 = 0.8602

celltype_label='clusters'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))
"""
Ductal: 0.9978
Ngn3 low EP: 0.9389
Ngn3 high EP: 0.9221
Beta: 0.8545
Alpha: 0.6965
Delta: 0.4857
Epsilon: 0.3099
"""

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4))

#################################################
#### velovi_wopreprocess
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import collections

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_velovi import *

method = 'velovi'
dataset_short = 'panINC'
dataset_long = 'pancreasINC'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/wopreprocess/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_v2.h5ad")
vae_split1 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_panINC_velovi_split1_wopreprocess_v2.pt', split1)

split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_v2.h5ad")
vae_split2 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_panINC_velovi_split2_wopreprocess_v2.pt', split2)

total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_v2.h5ad")
vae_total = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_panINC_velovi_total_wopreprocess_v2.pt', total)

raw = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
cell_index = np.array(np.where(raw.obs['clusters']!='Pre-endocrine')[0])
def create_adata_INC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

S_raw = raw.layers['spliced'][cell_index,:]
U_raw = raw.layers['unspliced'][cell_index,:]
raw = create_adata_INC(S=S_raw,U=U_raw,adata_old=raw)

total.uns['clusters_colors'] = split1.uns['clusters_colors'].copy()

## add velovi outputs to adata
add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)
## compute umap
compute_umap_pan(split1)
compute_umap_pan(split2)
compute_umap_pan(total)

split1.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_outputAdded_v2.h5ad")
split2.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_outputAdded_v2.h5ad")
total.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_outputAdded_v2.h5ad")

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

res = count_transition_matrix_pred(split1,split2,celltype_label='clusters')
np.sum(res), np.round(np.sum(res)/len(res),4)
# panINC+velovi_wopreprocess (2807, 0.9043)

celltype_label='clusters'
celltypes = split1.obs[celltype_label].cat.categories
celltypes_counter = collections.Counter(split1.obs[celltype_label])

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(ct+': '+str(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4)))
"""
Ductal: 0.9978
Ngn3 low EP: 0.9618
Ngn3 high EP: 0.9486
Beta: 0.9679
Alpha: 0.7401
Delta: 0.4286
Epsilon: 0.5211
"""

for ct in celltypes:
    idx = np.where(celltypes==ct)[0][0]
    print(np.round(np.sum(np.array(res)[split1.obs[celltype_label].values==celltypes[idx]])/celltypes_counter[celltypes[idx]], 4))



