import sctour as sct
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import test_timestep, print_message_with_time


split_seed = 317
method = 'sct'
dataset_long = 'larry'
dataset_short = 'larry'

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


times = []
means = []
medians = []
ptime_cors = []
latent_cors = []
for i in range(1,100):
    time = i/100
    times.append(time)
    mean_i,median_i,ptime_cor_i,latent_cor_i = test_timestep(adata_split1=split1,adata_split2=split2,adata_total=total,
                                                             tnode1=tnode_split1,tnode2=tnode_split2,tnode=tnode_total,time=time)
    means.append(mean_i)
    medians.append(median_i)
    ptime_cors.append(ptime_cor_i)
    latent_cors.append(latent_cor_i)
    print_message_with_time('################################ timestep='+str(i/100)+' done')

df = pd.DataFrame()
df['time'] = times
df['cos_sim_mean'] = means
df['cos_sim_median'] = medians
df['pseudotime_corr'] = ptime_cors
df['latenttime_corr'] = latent_cors
df.to_csv(data_folder+'larry_sct_velo_timestep.csv')

print_message_with_time('################################ Wrote data to '+data_folder+'larry_sct_velo_timestep.csv')


print_message_with_time('################################ All done')