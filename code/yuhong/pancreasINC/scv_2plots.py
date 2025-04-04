import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'scv'
dataset_long = 'pancreasINC'
dataset_short = 'panINC'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
raw = scv.read(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")

cell_index = np.array(np.where(raw.obs['clusters']!='Pre-endocrine')[0])
def create_adata_INC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}#adata_old.uns['clusters_colors'].copy()
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

S_raw = raw.layers['spliced'][cell_index,:]
U_raw = raw.layers['unspliced'][cell_index,:]
raw = create_adata_INC(S=S_raw,U=U_raw,adata_old=raw)

total.uns['clusters_colors'] = split1.uns['clusters_colors'].copy()

######################################################
## plot velocity
plot_velocity_scv_utv(adata_in=total,adata_raw=raw,fig_folder=fig_folder,fig_info="total",dataset=dataset_short,method=method)
plot_velocity_scv_utv(adata_in=split1,adata_raw=raw,fig_folder=fig_folder,fig_info="split1",dataset=dataset_short,method=method)
plot_velocity_scv_utv(adata_in=split2,adata_raw=raw,fig_folder=fig_folder,fig_info="split2",dataset=dataset_short,method=method)

plot_velocity_scv_utv(adata_in=total,adata_raw=raw,fig_folder=fig_folder,fig_info="recompF_total",dataset=dataset_short,method=method,recompute=False)
plot_velocity_scv_utv(adata_in=split1,adata_raw=raw,fig_folder=fig_folder,fig_info="recompF_split1",dataset=dataset_short,method=method,recompute=False)
plot_velocity_scv_utv(adata_in=split2,adata_raw=raw,fig_folder=fig_folder,fig_info="recompF_split2",dataset=dataset_short,method=method,recompute=False)

######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_hist_by_celltype(split1,split2,total,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## plot velo_conf
plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_veloConf_hist(adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## ptime
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total)
    scv.tl.velocity_pseudotime(split1)
    scv.tl.velocity_pseudotime(split2)

plot_pseudotime(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,method=method,fig_folder=fig_folder)

ptime_correlation_scatter_plot(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder)

# Spearman's corr
ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder,time_label='velocity_pseudotime')
ptime_correlation_scatter_spearman(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder,time_label='velocity_pseudotime')
ptime_correlation_scatter_spearman(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder,time_label='velocity_pseudotime')

# latent time
if not 'latent_time' in split1.obs.columns:
    scv.tl.latent_time(total)
    scv.tl.latent_time(split1)
    scv.tl.latent_time(split2)

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder,time_label='latent_time')
ptime_correlation_scatter_spearman(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder,time_label='latent_time')
ptime_correlation_scatter_spearman(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder,time_label='latent_time')

