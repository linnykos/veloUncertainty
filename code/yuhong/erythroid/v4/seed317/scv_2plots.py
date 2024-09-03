dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'scv'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

celltype_label = get_celltype_label(dataset_short)

total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_total_v4_outputAdded.h5ad')
split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split1_v4_outputAdded.h5ad')
split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split2_v4_outputAdded.h5ad')

## velocity
plot_velocity_scv_utv(adata_in=total,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label='celltype')
plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label='celltype')
plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label='celltype')

## cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)


c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
np.round(np.median(c1),4) # 0.5471
np.round(np.mean(c1),4) # 0.5238
np.round(np.var(c1),4) # 0.024

c2,n2 = compute_cosine_similarity_union(split1,split2,method)
np.round(np.median(c2),4) # 0.4315
np.round(np.mean(c2),4) # 0.4192
np.round(np.var(c2),4) # 0.0181

## confidence
plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)

np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) # 0.34453984

"""
>>> np.corrcoef([c2,total.obs['velocity_pseudotime'],total.obs['velocity_confidence']])
array([[ 1.        ,  0.29698886,  0.27629766],
       [ 0.29698886,  1.        , -0.33615143],
       [ 0.27629766, -0.33615143,  1.        ]])
"""
## pseudotime
plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label,ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed,celltype_label=celltype_label)


plot_pseudotime_diffusion(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime_diffusion(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
#plot_pseudotime_diffusion(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')

ptime_diffusion_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)

print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='velocity_pseudotime')
"""
Blood progenitors 1 0.4334
Blood progenitors 2 0.768
Erythroid1 0.2148
Erythroid2 0.2571
Erythroid3 0.0835
"""
print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='latent_time')


## latent time
plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)

# gene correlation
plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)

# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
np.round(np.mean(v2s_mean),4) # 0.2521
np.round(np.mean(v2s_median),4) # 0.2687

np.round(np.var(v2s_mean),4) # 
np.round(np.var(v2s_median),4) 

# paired: 0.4192 0.4315
# shuffled: 0.2521 0.2687

