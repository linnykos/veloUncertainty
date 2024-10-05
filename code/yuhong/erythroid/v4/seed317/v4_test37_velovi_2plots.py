import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_velovi import *

split_seed = 317
method = 'velovi'
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_test/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_test/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'


split1 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_37-3_v4.h5ad')
split2 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_37-7_v4.h5ad')

vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_37-3_v4.pt', split1)
vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_37-7_v4.pt', split2)

#######################################
## add velovi outputs to adata
print("############## Add velovi outputs to adata")

add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)

######################################################
## compute umap
print("############## Compute umap")

compute_umap_ery(split1)
compute_umap_ery(split2)

## write data
split1.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_37-3_v4_outputAdded.h5ad')
split2.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_37-7_v4_outputAdded.h5ad')

#split1 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_37-3_v4_outputAdded.h5ad')
#split2 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_37-7_v4_outputAdded.h5ad')
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)

#vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_37-3_v4.pt', split1)
#vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_37-7_v4.pt', split2)
vae_total = VELOVI.load('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_total_v4.pt', total)


######################################################
## plot velocity
plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="37-3",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="37-7",dataset=dataset_short,method=method,split_seed=split_seed)

######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) # 
c2,n2 = compute_cosine_similarity_union(split1,split2,method) # 
print(np.round(np.quantile(c1,[0.,.25,.5,.75,1.]),5))
print(np.round(np.quantile(c2,[0.,.25,.5,.75,1.]),5))


# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
print('shuffled mean and median')
print(np.round(np.mean(v2s_mean),4)) # 
print(np.round(np.mean(v2s_median),4)) # 

print(np.round(np.var(v2s_mean),4)) # 
print(np.round(np.var(v2s_median),4) )

# paired: 
# shuffled: 
