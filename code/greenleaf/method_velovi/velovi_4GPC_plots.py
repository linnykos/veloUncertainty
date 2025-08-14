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
method = 'velovi_GPC'
dataset_short = 'glf'
dataset_long = 'greenleaf'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'


total = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_total_GPC.h5ad') # 
split1 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split1_GPC.h5ad') # 
split2 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split2_GPC.h5ad') # 
vae_total = VELOVI.load(data_folder+'vae_glf_velovi_GPC_total_GPC.pt', total)
vae_split1 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split1_GPC.pt', split1)
vae_split2 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split2_GPC.pt', split2)

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

def compute_umap_GPC(adata):
    scv.pp.moments(adata, n_pcs=5, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=5) # used to be 10
    sc.tl.umap(adata)
    #scv.tl.velocity_graph(adata)

compute_umap_GPC(split1)
compute_umap_GPC(split2)
compute_umap(total, dataset_short) # n_neighbors=30, n_pcs=40

## write data
split1.write_h5ad(data_folder+'adata_glf_velovi_GPC_split1_GPC_outputAdded.h5ad')
split2.write_h5ad(data_folder+'adata_glf_velovi_GPC_split2_GPC_outputAdded.h5ad')
total.write_h5ad(data_folder+'adata_glf_velovi_GPC_total_GPC_outputAdded.h5ad')

"""
split1 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split1_GPC_outputAdded.h5ad')
split2 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split2_GPC_outputAdded.h5ad')
total = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_total_GPC_outputAdded.h5ad')

vae_total = VELOVI.load(data_folder+'vae_glf_velovi_GPC_total_GPC.pt', total)
vae_split1 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split1_GPC.pt', split1)
vae_split2 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split2_GPC.pt', split2)
"""

######################################################
## plot velocity
print('######## velocity')
plot_velocity(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
#plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
#plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)

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
    #scv.tl.velocity_confidence(split1)
    #scv.tl.velocity_confidence(split2)

#plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)

