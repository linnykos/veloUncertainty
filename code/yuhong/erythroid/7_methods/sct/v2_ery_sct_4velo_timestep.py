import sctour as sct
import scvelo as scv
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
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split2_v2.h5ad') # 
tnode_total = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_total_v2.pth')
tnode_split1 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_seed317_split1_v2.pth')
tnode_split2 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_seed317_split2_v2.pth')
print('######### Read data done!')

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
    return np.round(np.mean(cos_sim),6),np.round(np.median(cos_sim),6),ptime_cor[0,1]

times = []
means = []
medians = []
ptime_cors = []
for i in range(1,91):
    time = i/100
    times.append(time)
    mean_i,median_i,cor_i = test_timestep(adata_split1=split1,adata_split2=split2,adata_total=total,
                                          tnode1=tnode_split1,tnode2=tnode_split2,tnode=tnode_total,time=time)
    means.append(mean_i)
    medians.append(median_i)
    ptime_cors.append(cor_i)
    print(str(time)+' done')

df = pd.DataFrame()
df['time'] = times
df['cos_sim_mean'] = means
df['cos_sim_median'] = medians
df['pseudotime_corr'] = ptime_cors
df.to_csv(data_folder+'v2_'+dataset_long+'/'+method+'/'+dataset_short+'_'+method+'velo_timestep.csv')
print('######### Wrote data to '+data_folder+'v2_'+dataset_long+'/'+method+'/'+dataset_short+'_'+method+'velo_timestep.csv')