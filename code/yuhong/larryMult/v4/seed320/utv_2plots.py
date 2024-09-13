dataset_long = 'larryMult'
dataset_short = 'larryMult'
method = 'utv'
split_seed=320

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

total = read_data_v4(dataset_long,dataset_short,method,split_seed=317,data_version='total',allgenes=False,outputAdded=True)

split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)

colors = ["#6e8ea1","#ffab6e","#dba8bc","#a0a0a0","#c4c88a","#87c3c9"]
split1.uns['state_info_colors'] = colors
split2.uns['state_info_colors'] = colors

compute_umap(split1, dataset_short)
compute_umap(split2, dataset_short)

split1.obsm['X_umapOriginal'] = split1.obsm['X_umap'].copy()
split1.obsm['X_umapOriginal'][:,0] = np.array(split1.obs['SPRING-x'])
split1.obsm['X_umapOriginal'][:,1] = np.array(split1.obs['SPRING-y'])

split2.obsm['X_umapOriginal'] = split2.obsm['X_umap'].copy()
split2.obsm['X_umapOriginal'][:,0] = np.array(split2.obs['SPRING-x'])
split2.obsm['X_umapOriginal'][:,1] = np.array(split2.obs['SPRING-y'])

#split1.write_h5ad(data_folder+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
#split2.write_h5ad(data_folder+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
#total.write_h5ad(data_folder+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad')
"""
total = read_data_v4(dataset_long,dataset_short,method,data_version='total',allgenes=False,outputAdded=True)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
"""
## velocity
plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed)

## cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) # 868
c2,n2 = compute_cosine_similarity_union(split1,split2,method) # 1367
np.quantile(c1,[0.,.25,.5,.75,1.]) # [-0.80632752,  0.70185971,  0.85641897,  0.93676168,  0.96146846]
np.quantile(c2,[0.,.25,.5,.75,1.]) # [-0.78773004,  0.68730429,  0.83893331,  0.91563132,  0.93775107]

######################################################
## plot velo_conf
if (not 'velocity_confidence' in split1.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)

np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) # 0.04445349

######################################################
## ptime
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)

plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed,celltype_label='state_info')
# 0.788

if not 'latent_time' in split1.obs.columns:
    scv.tl.latent_time(total)
    scv.tl.latent_time(split1)
    scv.tl.latent_time(split2)

plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
np.round(np.mean(v2s_mean),4) # 
np.round(np.mean(v2s_median),4) # 

np.round(np.var(v2s_mean),4) # 
np.round(np.var(v2s_median),4) # 

# paired: 0.7049 0.8389
# shuffled: 0.1601 0.3044

np.quantile([np.var(split1.layers['velocity'][:,i]) for i in range(split1.shape[1])], [0.,.25,.5,.75,1.])
# array([7.91685125e-06, 1.11124450e-02, 2.83529945e-02, 1.90240689e-01,
#        2.15817334e+03])
np.quantile([np.var(split2.layers['velocity'][:,i]) for i in range(split2.shape[1])], [0.,.25,.5,.75,1.])
# array([1.30304642e-07, 7.15183408e-03, 1.93593809e-02, 7.70638399e-02,
#        2.41481714e+03])
