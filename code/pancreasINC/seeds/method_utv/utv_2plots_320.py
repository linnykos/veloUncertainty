dataset_long = 'pancreasINC'
dataset_short = 'panINC'
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

total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)

compute_umap(total,dataset_short)
compute_umap(split1,dataset_short)
compute_umap(split2,dataset_short)

raw = read_raw_adata(dataset_short)
split1.uns['clusters_colors'] = raw.uns['clusters_colors'].copy()
split2.uns['clusters_colors'] = raw.uns['clusters_colors'].copy()
total.uns['clusters_colors'] = raw.uns['clusters_colors'].copy()

"""

#total.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad')
#split1.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
#split2.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')


#total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
#split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
#split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
"""

## velocity
#plot_velocity_scv_utv(adata_in=total,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed)

## cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)


c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
# 1497, 2503

np.quantile(c1, [0.,.25,.5,.75,1.]) # [-0.04854831,  0.70391344,  0.811932  ,  0.9531423 ,  0.97302145]
np.quantile(c2, [0.,.25,.5,.75,1.]) # [-0.04832041,  0.70169404,  0.80931835,  0.9450743 ,  0.96540783]


## confidence
if (not 'velocity_confidence' in split1.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)

np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) # 0.19600904


######################################################
## ptime
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)

plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
#plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)
# 0.877

if not 'latent_time' in split1.obs.columns:
    scv.tl.recover_dynamics(total,n_jobs=8)
    scv.tl.velocity_graph(total,n_jobs=8)
    scv.tl.latent_time(total)
    scv.tl.recover_dynamics(split1,n_jobs=8)
    scv.tl.velocity_graph(split1,n_jobs=8)
    scv.tl.latent_time(split1)
    scv.tl.recover_dynamics(split2,n_jobs=8)
    scv.tl.velocity_graph(split2,n_jobs=8)
    scv.tl.latent_time(split2)

plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)
# 0.89

np.corrcoef([c2,total.obs['velocity_confidence'],total.obs['velocity_pseudotime'],total.obs['latent_time']])
"""
array([[ 1.        , -0.23104282, -0.19566166, -0.58183944],
       [-0.23104282,  1.        ,  0.47805072,  0.54313432],
       [-0.19566166,  0.47805072,  1.        ,  0.77345141],
       [-0.58183944,  0.54313432,  0.77345141,  1.        ]])
"""

# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
( np.round(np.mean(v2s_mean),4) , np.round(np.mean(v2s_median),4) )

np.round(np.var(v2s_mean),4) # 
np.round(np.var(v2s_median),4) 
np.round(np.var(v2s_mean),4) # 
np.round(np.var(v2s_median),4) # 

# paired: 0.7841 0.8093
# shuffled: 0.2679, 0.4556

######################################
## random, cosine similarity by cells and neighborhoods

total.obs['clusters']

#
import collections
celltype_label = 'clusters'
idx_mix = []
idx_mid = []
for i in range(total.shape[0]):
    idx = np.where(total.obsp['connectivities'][i].todense()>0)[1]
    celltype_i = total.obs[celltype_label][i]
    celltype_nb = total.obs[celltype_label][idx]
    counter = collections.Counter(total.obs[celltype_label][idx])
    if len(counter)==1 and counter[celltype_i]>0:
        idx_mid.append(i)
    else:
        idx_mix.append(i)

[ len(idx_mid),len(idx_mix) ] # [1239, 1865]
# cell index
idx_low_cos = np.where(c2 < np.quantile(c2,[.25])[0])[0]

len(np.intersect1d(idx_low_cos,idx_mid)) # 559
len(np.intersect1d(idx_low_cos,idx_mix)) # 217

np.min(c2[idx_low_cos])
np.min(c2[idx_mid])
np.min(c2[idx_mix])

collections.Counter(total.obs[celltype_label][np.where(c2 < np.quantile(c2,[.25])[0])[0]])
# Counter({'Ngn3 high EP': 328, 'Beta': 313, 'Alpha': 94, 'Epsilon': 27, 'Delta': 14})

i = 0
j = np.where(c2 < np.quantile(c2,[.25])[0])[0][i]
total.obs[celltype_label][j]
collections.Counter(total.obs[celltype_label][np.where(total.obsp['connectivities'][j].todense()>0)[0]])

for i in np.where(c2 < np.quantile(c2,[.25])[0])[0]:
    counter_nb = collections.Counter(total.obs[celltype_label][np.where(total.obsp['connectivities'][i].todense()>0)[0]])
    print(total.obs[celltype_label][i],counter_nb[total.obs[celltype_label][i]],len(counter_nb))

for i in np.where(c2 > np.quantile(c2,[.9])[0])[0]:
    counter_nb = collections.Counter(total.obs[celltype_label][np.where(total.obsp['connectivities'][i].todense()>0)[0]])
    print(total.obs[celltype_label][i],counter_nb[total.obs[celltype_label][i]],len(counter_nb))

for i in np.where(c2 > np.quantile(c2,[.25])[0])[0]:
    if c2[i] > np.quantile(c2,[.5])[0]: continue
    else:
        counter_nb = collections.Counter(total.obs[celltype_label][np.where(total.obsp['connectivities'][i].todense()>0)[0]])
        print(total.obs[celltype_label][i],
            counter_nb[total.obs[celltype_label][i]],len(counter_nb))

np.quantile(c2[idx_mid],[0.,.25,.5,.75,1.])
np.quantile(c2[idx_mix],[0.,.25,.5,.75,1.])
np.quantile(c2,[0.,.25,.5,.75,1.])
"""
"""