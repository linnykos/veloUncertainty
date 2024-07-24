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
dataset_long = 'pancreasINC'
dataset_short = 'panINC'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 

tnode_split1 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split1_v2.pth')
tnode_split2 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split2_v2.pth')
print_message_with_time('######### Read data done!')

def test_timestep_split(adata_in,tnode,time):
    adata = adata_in.copy()
    adata.layers['velocity'] = compute_sctour_velocity(tnode, timestep=time)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    ptime_cor = np.corrcoef(adata.obs['ptime'],adata.obs['velocity_pseudotime'])
    print(ptime_cor[0,1])
    return ptime_cor[0,1]

times = []
ptime_cor1 = []
ptime_cor2 = []
for i in range(1,1000):
    time = i/1000
    times.append(time)
    cor1 = test_timestep_split(adata_in=split1,tnode=tnode_split1,time=time)
    cor2 = test_timestep_split(adata_in=split2,tnode=tnode_split2,time=time)
    print('*** ptime split1: ')
    ptime_cor1.append(cor1)
    print('*** ptime split2: ')
    ptime_cor2.append(cor2)
    print_message_with_time('######### '+str(i)+' done')

test_timestep_split(adata_in=split1,tnode=tnode_split1,time=1/20)
test_timestep_split(adata_in=split2,tnode=tnode_split2,time=1/20)

df = pd.DataFrame()
df['time'] = times
df['ptime_cor1'] = ptime_cor1
df['ptime_cor2'] = ptime_cor2
df.to_csv(data_folder+'v2_'+dataset_long+'/'+method+'/'+dataset_short+'_'+method+'_velo_timestep_splits.csv')
print_message_with_time('######### Wrote data to '+data_folder+'v2_'+dataset_long+'/'+method+'/'+dataset_short+'_'+method+'velo_timestep_splits.csv')