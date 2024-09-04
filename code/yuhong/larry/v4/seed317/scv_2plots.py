dataset_long = "larry"
dataset_short = "larry"
method = "scv"
split_seed = 317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
adata_prefix = 'adata_'+dataset_short+'_'+method

fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4.h5ad')
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4.h5ad')
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4.h5ad')

raw = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/larry.h5ad') # n_obs × n_vars = 49302 × 23420
raw.obsm['X_umap']=raw.obsm['X_emb']

colors = ["#6e8ea1","#ffab6e","#91b68a","#d18887","#a892c4","#b2967e","#dba8bc","#a0a0a0","#c4c88a","#87c3c9","#ffcccc"]
split1.uns['state_info_colors'] = colors
split2.uns['state_info_colors'] = colors
total.uns['state_info_colors'] = colors

split1.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
split2.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
total.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()

######################################################
## plot velocity
plot_velocity_scv_utv(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)

######################################################
## plot cosine similarity
#c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
#c2,n2 = compute_cosine_similarity_union(split1,split2,method)
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
#plot_cosine_similarity_violinplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

######################################################
## plot velo_conf
if (not method=='sct') and (not 'velocity_confidence' in total.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,raw,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)

np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) # 0.29408268

######################################################
## ptime
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total)
    scv.tl.velocity_pseudotime(split1)
    scv.tl.velocity_pseudotime(split2)

plot_pseudotime(adata_in=split1,fig_name='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,fig_name='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=total,fig_name='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed,celltype_label='state_info')

if not 'latent_time' in split1.obs.columns:
    scv.tl.latent_time(total)
    scv.tl.latent_time(split1)
    scv.tl.latent_time(split2)

plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')