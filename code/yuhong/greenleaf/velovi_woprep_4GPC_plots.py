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
method = 'velovi_woprep_GPC'
dataset_short = 'glf'
dataset_long = 'greenleaf'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'


total = sc.read_h5ad(data_folder+'adata_glf_'+method+'_total_GPC.h5ad') # 
split1 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split1_GPC.h5ad') # 
split2 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split2_GPC.h5ad') # 
vae_total = VELOVI.load(data_folder+'vae_glf_'+method+'_total_GPC.pt', total)
vae_split1 = VELOVI.load(data_folder+'vae_glf_'+method+'_split1_GPC.pt', split1)
vae_split2 = VELOVI.load(data_folder+'vae_glf_'+method+'_split2_GPC.pt', split2)

total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()

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
split1.write_h5ad(data_folder+'adata_glf_'+method+'_split1_GPC_outputAdded.h5ad')
split2.write_h5ad(data_folder+'adata_glf_'+method+'_split2_GPC_outputAdded.h5ad')
total.write_h5ad(data_folder+'adata_glf_'+method+'_total_GPC_outputAdded.h5ad')

"""
split1 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split1_GPC_outputAdded.h5ad')
split2 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split2_GPC_outputAdded.h5ad')
total = sc.read_h5ad(data_folder+'adata_glf_'+method+'_total_GPC_outputAdded.h5ad')

vae_total = VELOVI.load(data_folder+'vae_glf_'+method+'_total_GPC.pt', total)
vae_split1 = VELOVI.load(data_folder+'vae_glf_'+method+'_split1_GPC.pt', split1)
vae_split2 = VELOVI.load(data_folder+'vae_glf_'+method+'_split2_GPC.pt', split2)
"""

######################################################
## plot velocity
print('######## velocity')
plot_velocity(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)

######################################################
## plot cosine similarity
print('######## cosine similarity')
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) 
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
np.quantile(c1,[0.,.25,.5,.75,1.]) 
np.quantile(c2,[0.,.25,.5,.75,1.]) 

######################################################
## plot velo_conf
print('######## velocity confidence')
if (not 'velocity_confidence' in total.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)

print('correlation of confidence between splits')
np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) 

######################################################
## ptime
print('######## pseudotime')
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)

plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)







