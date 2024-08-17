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

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from sctour_misc import *
from v2_functions import *

method = 'sct'
dataset_long = 'larry'
dataset_short = 'larry'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
tnode_split1 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split1_v2.pth')

print_message_with_time('################################ Read data done')

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
for i in range(1,100):
    time = i/100
    times.append(time)
    print('*** ptime split1: ')
    cor1 = test_timestep_split(adata_in=split1,tnode=tnode_split1,time=time)
    ptime_cor1.append(cor1)
    print_message_with_time('######### '+str(i)+' done')

df = pd.DataFrame()
df['time'] = times
df['ptime_cor1'] = ptime_cor1
df.to_csv(data_folder+'v2_'+dataset_long+'/'+method+'/'+dataset_short+'_'+method+'_velo_timestep_split1.csv')
print_message_with_time('######### Wrote data to '+data_folder+'v2_'+dataset_long+'/'+method+'/'+dataset_short+'_'+method+'velo_timestep_split1.csv')
