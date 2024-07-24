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
from v2_functions import *
from sctour_misc import *

method = 'sct'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split2_v2.h5ad') # 
#raw = sc.read_h5ad
#raw.obsm['X_umap']=raw.obsm['X_emb']
tnode_total = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_total_v2.pth')
tnode_split1 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_seed317_split1_v2.pth')
tnode_split2 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_seed317_split2_v2.pth')

v_1_100 = compute_sctour_velocity(tnode_total, timestep=1/100)
v1_1_100 = compute_sctour_velocity(tnode_split1, timestep=1/100)
v2_1_100 = compute_sctour_velocity(tnode_split2, timestep=1/100)

cos_sim_1_100=np.diag(cosine_similarity(v1_1_100,v1_1_100))

v1_1_100-split1.layers['velocity']

np.diag(cosine_similarity(split1.layers['velocity'],split2.layers['velocity']))


cos_sim_1_100,Ngenes_1_100 = compute_cosine_similarity(split1,split2,method='sct')

np.quantile(cos_sim_1_10,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
np.quantile(cos_sim_1_100,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

"""
np.quantile(cos_sim_1_10,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
array([-0.80406284,  0.86307514,  0.91327848,  0.93471146,  0.94986193,
        0.96042544,  0.96817834,  0.97405195,  0.97895435,  0.98359616,
        0.99198878])
np.quantile(cos_sim_1_100,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
array([-0.85378456,  0.81818457,  0.88342689,  0.91475611,  0.93092065,
        0.94250852,  0.95063586,  0.95685688,  0.96255744,  0.9682582 ,
        0.98711741])
"""
import scvelo as scv
# time=1/100
total.obs['ptime']
scv.tl.velocity_graph(total)
scv.tl.velocity_pseudotime(total)
total.obs['velocity_pseudotime']
np.corrcoef(total.obs['ptime'],total.obs['velocity_pseudotime']) # 0.75022968

# time=1/10
## total
total_1_10 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
v_1_10 = compute_sctour_velocity(tnode_total, timestep=1/10)
total_1_10.layers['velocity'] = v_1_10
scv.tl.velocity_graph(total_1_10)
scv.tl.velocity_pseudotime(total_1_10)
np.corrcoef(total_1_10.obs['ptime'],total_1_10.obs['velocity_pseudotime']) # 0.8012634
## splits
split1_1_10 = split1.copy()
del split1_1_10.layers['velocity']
split2_1_10 = split2.copy()
del split2_1_10.layers['velocity']
split1_1_10.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=1/10)
split2_1_10.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=1/10)
cos_sim_1_10,Ngenes_1_10 = compute_cosine_similarity(split1_1_10,split2_1_10,method='sct')

def test_timestep_split(adata_in,tnode,time):
    adata = adata_in.copy()
    adata.layers['velocity'] = compute_sctour_velocity(tnode, timestep=time)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    ptime_cor = np.corrcoef(adata.obs['ptime'],adata.obs['velocity_pseudotime'])
    print(ptime_cor[0,1])

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
    return cos_sim,ptime_cor

# time=1/2
cos_sim_1_2,ptime_cor_1_2 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/2)
np.quantile(cos_sim_1_2,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.])
[np.mean(cos_sim_1_2),np.median(cos_sim_1_2)] # [0.900714, 0.93990314]
ptime_cor_1_2 # 0.58324394
"""
[0.900714, 0.93990314]
0.58324394
"""
test_timestep_split(adata_in=split1,tnode=tnode_split1,time=1/2) # 0.8225269723618918
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=1/2) # 0.8084339958720952

# time=5/12 (=10/24)
test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=5/12)
"""
[0.91310376, 0.9518938]
0.48160690069151296
"""
# time=9/24
test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=9/24)
"""
[0.9178761, 0.9566163]
0.41912900420270377
"""
# time=1/3 (=4/12=8/24)
test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/3)
"""
[0.92153114, 0.96040744]
0.3946008691075687
"""
test_timestep_split(adata_in=split1,tnode=tnode_split1,time=1/3) # 0.8213109619120375
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=1/3) # 0.8118910956637128

# time=7/24 (=14/48)
test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=7/24)
"""
[0.92394394, 0.9634903]
0.5037084533564883
"""
# time=13/48 
test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=13/48)
"""
[0.9246232, 0.9644784]
0.6161328799620464
"""

# time=1/4 (=3/12=6/24)
test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/4)
"""
[0.9249223, 0.9652108]
0.7071871787435123
"""
test_timestep_split(adata_in=split1,tnode=tnode_split1,time=1/4) # 0.8209524300049058
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=1/4) # 0.8136089214798514

# time=9/40
test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=9/40)

test_timestep_split(adata_in=split1,tnode=tnode_split1,time=9/40) 
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=9/40)

# time=1/5
cos_sim_1_5,ptime_cor_1_5 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/5)
print([np.mean(cos_sim_1_5), np.median(cos_sim_1_5)]) # [0.92421585, 0.96555763]
print(ptime_cor_1_5[0,1]) # 0.7792453538490794
"""
[0.92421585, 0.96555763] 
0.7792453538490794
"""
test_timestep_split(adata_in=split1,tnode=tnode_split1,time=1/5) # 0.821367110294137
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=1/5) # 0.8159337576625478

# time=1/10
cos_sim_1_10,ptime_cor_1_10 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/10)
"""
[0.91560787, 0.96042544]
0.7684178268965088
"""
test_timestep_split(adata_in=split1,tnode=tnode_split1,time=1/10) # 0.824713879312453
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=1/10) # 0.8213970375026902

# time=1/20
cos_sim_1_20,ptime_cor_1_20 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/20)
"""
[0.90546733, 0.95265067]
0.7576828946960318
"""
test_timestep_split(adata_in=split1,tnode=tnode_split1,time=1/20) # 0.8280501058711189
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=1/20) # 0.8249701683240939

# time=1/50
cos_sim_1_50,ptime_cor_1_50 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/50)
"""
[0.8974141, 0.9452971]
0.7502663703710624
"""

# time=1/100
cos_sim_1_100,ptime_cor_1_100 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/100)
"""
[0.89437264, 0.9425085]
0.7502371783928625
"""

# time=1/200
cos_sim_1_200,ptime_cor_1_200 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/200)
"""
[0.8927245, 0.9406786]
0.7460500148492272
"""
# time=1/250
cos_sim_1_250,ptime_cor_1_250 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/250)
"""
[0.8924052, 0.9402451]
0.7450498698173291
"""
# time=1/500
cos_sim_1_500,ptime_cor_1_500 = test_timestep(split1,split2,total,tnode_split1,tnode_split2,tnode_total,time=1/500)
"""
[0.8917709, 0.93950003]
0.7422150507668348
"""
