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
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

"""
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

compute_umap_ery(total)
compute_umap_ery(split1)
compute_umap_ery(split2)

compute_umap(split1, dataset_short)
compute_umap(split2, dataset_short)
compute_umap(total, dataset_short)

## write data
split1.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
split2.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
total.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad')

"""
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

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) # 1480
c2,n2 = compute_cosine_similarity_union(split1,split2,method) # 2520
np.quantile(c1,[0.,.25,.5,.75,1.]) # [0.63100588, 0.91638279, 0.94828916, 0.96809444, 0.99312139]
np.quantile(c2,[0.,.25,.5,.75,1.]) # [0.62972143, 0.91227215, 0.94539505, 0.96535868, 0.99224501]

######################################################
## plot velo_conf
if (not 'velocity_confidence' in total.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)

np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) # 0.43612837

######################################################
## ptime
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)

plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)
# 0.839

"""
if not 'latent_time' in split1.obs.columns:
    scv.tl.latent_time(total)
    scv.tl.latent_time(split1)
    scv.tl.latent_time(split2)

# KeyError: 'fit_likelihood'
plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
"""

# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
np.round(np.mean(v2s_mean),4) # 0.8757
np.round(np.mean(v2s_median),4) # 0.8994

np.round(np.var(v2s_mean),4) # 
np.round(np.var(v2s_median),4) 

# paired: 0.9313 0.9454
# shuffled: 0.8757 0.8994

np.quantile([np.var(total.layers['velocity'][:,i]) for i in range(2000)],[0.,.25,.5,.75,1.])
np.quantile([np.var(split1.layers['velocity'][:,i]) for i in range(2000)],[0.,.25,.5,.75,1.])
np.quantile([np.var(split2.layers['velocity'][:,i]) for i in range(2000)],[0.,.25,.5,.75,1.])

# [0.00000000e+00, 7.64297667e-07, 1.44990809e-05, 3.13360120e-04, 8.36765747e+01]
# [0.00000000e+00, 3.34343170e-07, 3.48852461e-06, 3.42113917e-05, 3.72625694e+01]
# [0.00000000e+00, 4.86066682e-07, 4.30251725e-06, 3.82065136e-05, 3.74937248e+01]
