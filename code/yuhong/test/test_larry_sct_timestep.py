import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity

import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from sctour_misc import *
from v2_functions import *

method = 'sct'
dataset_long = 'larry'
dataset_short = 'larry'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
"""
raw = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/larry.h5ad') # n_obs × n_vars = 49302 × 23420
raw.obsm['X_umap']=raw.obsm['X_emb']

import matplotlib as mpl
def rgb2hex(rgb):
    r = int(rgb[0]*255)
    g = int(rgb[1]*255)
    b = int(rgb[2]*255)
    return "#{:02x}{:02x}{:02x}".format(r,g,b)

split1.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]
split2.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]
total.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]
"""
tnode_total = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_total_v2.pth')

tnode_split1 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split1_v2.pth')
tnode_split2 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split2_v2.pth')

velo_split1 = compute_sctour_velocity(tnode_split1, timestep=1/10)
velo_split2 = compute_sctour_velocity(tnode_split2, timestep=1/10)

import scvelo as scv
def test_timestep(adata_split1,adata_split2,adata_total,tnode1,tnode2,tnode,time):
    split1 = adata_split1.copy()
    split2 = adata_split2.copy()
    total = adata_total.copy()
    split1.layers['velocity'] = compute_sctour_velocity(tnode1, timestep=time)
    split2.layers['velocity'] = compute_sctour_velocity(tnode2, timestep=time)
    total.layers['velocity'] = compute_sctour_velocity(tnode, timestep=time)
    cos_sim,Ngenes = compute_cosine_similarity(split1,split2,method='sct')
    scv.tl.velocity_graph(total)
    scv.tl.velocity_pseudotime(total)
    ptime_cor = np.corrcoef(total.obs['ptime'],total.obs['velocity_pseudotime'])
    print([np.mean(cos_sim), np.median(cos_sim)])
    print(ptime_cor[0,1])


np.mean(compute_cosine_similarity(split1,split2,method='sct')[0]) # 0.6707662
np.median(compute_cosine_similarity(split1,split2,method='sct')[0]) # 0.81822276

test_timestep(adata_split1=split1,adata_split2=split2,adata_total=total,tnode1=tnode_split1,tnode2=tnode_split2,tnode=tnode_total,time=1/100)
"""
[0.6707662, 0.81822276]
0.324890113967215
"""

for time in [1/2,1/3,1/4,1/5,1/10,1/25,1/50,1/100,1/200,1/500]: 
    print('###################################################################### ')
    print('################################### '+str(np.round(time,4)))
    test_timestep(adata_split1=split1,adata_split2=split2,adata_total=total,tnode1=tnode_split1,tnode2=tnode_split2,tnode=tnode_total,time=time)
    print('Time='+str(np.round(time,4))+' done!')


"""
time=1/2
[0.58148944, 0.6294818]
0.5204163967537302
(array([0.36991194, 0.6122273 , 0.59724903, ..., 0.9173325 , 0.86590976, 0.55407137], dtype=float32), 
 array([[1.       , 0.5204164],
       [0.5204164, 1.       ]]))

time=1/3
[0.6474253, 0.74257505]
0.5852286822170534
(array([0.8085468 , 0.67628473, 0.2146571 , ..., 0.93948615, 0.72326475, 0.5986347 ], dtype=float32), 
 array([[1.        , 0.58522868],
       [0.58522868, 1.        ]]))

time=1/4
[0.67303485, 0.79298115]
0.5923016506166456
(array([0.8414159 , 0.6570008 , 0.12108883, ..., 0.9450992 , 0.63772404, 0.6715509 ], dtype=float32), 
 array([[1.        , 0.59230165],
       [0.59230165, 1.        ]]))

time=1/5
[0.686955, 0.8243335]
0.5966680118787147
(array([0.8878193 , 0.62751204, 0.11528225, ..., 0.93498105, 0.62352127, 0.6367445 ], dtype=float32), 
 array([[1.        , 0.59666801],
       [0.59666801, 1.        ]]))

time=1/10
[0.6758771, 0.8325791]
0.315832026189503
(array([0.8173126 , 0.46920297, 0.76714563, ..., 0.40456256, 0.75751305, 0.9243948 ], dtype=float32), 
 array([[1.        , 0.31583203],
       [0.31583203, 1.        ]]))

time=1/25
[0.668259, 0.8243551]
0.31239882321774054
(array([0.9258586 , 0.28131756, 0.8998147 , ..., 0.78898704, 0.8471142 , 0.89710903], dtype=float32), 
 array([[1.        , 0.31239882],
       [0.31239882, 1.        ]]))

time=1/50
[0.6709529, 0.8205998]
0.3296690655616754
(array([0.9341025 , 0.16277966, 0.90364444, ..., 0.7596393 , 0.8582997 , 0.8794918 ], dtype=float32), 
 array([[1.        , 0.32966907],
       [0.32966907, 1.        ]]))

time=1/100
[0.6707662, 0.81822276]
0.3248869451883729
(array([0.91538584, 0.1110096 , 0.90529454, ..., 0.7863237 , 0.86207354, 0.86849666], dtype=float32), 
 array([[1.        , 0.32488695],
       [0.32488695, 1.        ]]))

time=1/200
[0.67021155, 0.8171601]
0.32420737645377384
(array([0.89629304, 0.10935047, 0.9062141 , ..., 0.7820419 , 0.86047584, 0.8688125 ], dtype=float32), 
 array([[1.        , 0.32420738],
       [0.32420738, 1.        ]]))

time=1/500
[0.6699041, 0.81655073]
0.3247905438219868
(array([0.8961598 , 0.08163041, 0.9060961 , ..., 0.78087187, 0.84762967, 0.86801165], dtype=float32), 
 array([[1.        , 0.32479054],
       [0.32479054, 1.        ]]))
"""

split1.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=1/5)
split2.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=1/5)
scv.tl.velocity_graph(split1)
scv.tl.velocity_pseudotime(split1)
scv.tl.velocity_graph(split2)
scv.tl.velocity_pseudotime(split2)
np.corrcoef(split1.obs['velocity_pseudotime'],split2.obs['velocity_pseudotime']) # 0.0951427
np.corrcoef(split1.obs['ptime'],split2.obs['ptime']) # 0.92236945
cos_sim, Ngenes=compute_cosine_similarity(split1,split2,method='sct') 
np.quantile(cos_sim,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.])
# array([-0.87171334,  0.20372799,  0.49924864,  0.65870616,  0.76200461, 
#         0.82433349,  0.86604184,  0.89874809,  0.92549888,  0.95067748, 0.99357533])
np.mean(cos_sim) # 0.686955
np.median(cos_sim) # 0.8243335
