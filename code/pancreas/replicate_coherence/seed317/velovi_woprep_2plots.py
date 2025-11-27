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
method = 'velovi_woprep'
dataset_short = 'pan'
dataset_long = 'pancreas'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split1_v4.pt', split1)

split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split2_v4.pt', split2)

total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
vae_total = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_total_v4.pt', total)

#######################################
## add velovi outputs to adata
print_message_with_time("############## Add velovi outputs to adata")

add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)

######################################################
## compute umap
print_message_with_time("############## Compute umap")

compute_umap(split1, dataset_short)
compute_umap(split2, dataset_short)
compute_umap(total, dataset_short)

## write data
split1.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
split2.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
total.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad')


### read the outputAdded version
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)

vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split1_v4.pt', split1)
vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split2_v4.pt', split2)
vae_total = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_total_v4.pt', total)

######################################################
# gene correlation
plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)

######################################################
## plot velocity
plot_velocity(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)

######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) # 1487
c2,n2 = compute_cosine_similarity_union(split1,split2,method) # 2513
np.quantile(c1,[0.,.25,.5,.75,1.]) # [0.2841692 , 0.69088823, 0.74792105, 0.79119323, 0.94226491]
np.quantile(c2,[0.,.25,.5,.75,1.]) # [0.27720253, 0.67109071, 0.72246853, 0.76629862, 0.93210125]

######################################################
## plot velo_conf
if (not 'velocity_confidence' in total.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)

######################################################
# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
np.round(np.mean(v2s_mean),4) # 0.4052
np.round(np.mean(v2s_median),4) # 0.4139

np.round(np.var(v2s_mean),4) # 
np.round(np.var(v2s_median),4) 

# paired: 0.7122 0.7225
# shuffled: 0.4052 0.4139

