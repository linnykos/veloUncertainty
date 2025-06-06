split_seed = 317
method = 'sct'
dataset_long = 'pancreas'
dataset_short = 'pan'

import sctour as sct
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import print_message_with_time
from sctour_misc import *


import torch
import random
from scipy.stats import spearmanr
def test_timestep(adata_split1,adata_split2,adata_total,tnode1,tnode2,tnode,time,sct_seed=615):
    total = adata_total.copy()
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    total.layers['velocity'] = compute_sctour_velocity(tnode, timestep=time)
    scv.tl.velocity_graph(total,n_jobs=8)
    scv.tl.velocity_pseudotime(total)
    ptime_cor = spearmanr(total.obs['ptime'], total.obs['velocity_pseudotime']).correlation
    print(ptime_cor)
    return np.round(ptime_cor,5)


data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4.h5ad') # 
tnode_total = sct.predict.load_model(data_folder+tnode_prefix+'_total_v4.pth')
tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')

print_message_with_time('################################ Read data done')

scv.pp.moments(total, n_pcs=30, n_neighbors=30)
sc.tl.pca(total, svd_solver="arpack")
sc.pp.neighbors(total, n_neighbors=30, n_pcs=40)

times = []
ptime_cors = []
for i in range(1,100):
    time = i/100
    times.append(time)
    ptime_cor_i = test_timestep(adata_split1=split1,adata_split2=split2,adata_total=total,
                                tnode1=tnode_split1,tnode2=tnode_split2,tnode=tnode_total,time=time)

    ptime_cors.append(ptime_cor_i)
    print_message_with_time('################################ timestep='+str(i/100)+' done')

df = pd.DataFrame()
df['time'] = times
df['pseudotime_corr'] = ptime_cors
df.to_csv(data_folder+dataset_short+'_sct_velo_timestep.csv')

print_message_with_time('################################ Wrote data to '+data_folder+dataset_short+'_sct_velo_timestep.csv')

print_message_with_time('################################ All done')
print(df.iloc[np.where(df['pseudotime_corr']==np.max(df['pseudotime_corr']))[0][0]]['time'])

"""
#########################################
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
df = pd.read_csv(data_folder+dataset_short+'_sct_velo_timestep.csv')
df.iloc[np.where(df['pseudotime_corr']==np.max(df['pseudotime_corr']))[0][0]]['time'] 
"""