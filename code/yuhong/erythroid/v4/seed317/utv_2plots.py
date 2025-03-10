dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'utv'
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
from v4_functions_scv import plot_velocity_scv_utv
from v4_functions import *

celltype_label = get_celltype_label(dataset_short)
"""
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)

compute_umap_ery(total)
compute_umap_ery(split1)
compute_umap_ery(split2)

total.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad')
split1.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
split2.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
"""

total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)

## velocity
plot_velocity_scv_utv(adata_in=total,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)

## cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)


c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)

np.round(np.median(c1),4) # 0.9967
np.round(np.mean(c1),4) # 0.8613
np.round(np.var(c1),4) # 0.1961

np.round(np.median(c2),4) # 0.9949
np.round(np.mean(c2),4) # 0.8597
np.round(np.var(c2),4) # 0.1953

np.quantile(c1, [0.,.25,.5,.75,1.]) # [-0.99336785,  0.99321961,  0.9967162 ,  0.99697784,  0.99711633]
np.quantile(c2, [0.,.25,.5,.75,1.]) # [-0.99159464,  0.99136127,  0.99493376,  0.99520167,  0.99534509]


## confidence
if (not 'velocity_confidence' in split1.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=celltype_label)

np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) # 0.29686059


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
# 

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
# 0.886

np.corrcoef([c2,total.obs['velocity_confidence'],total.obs['velocity_pseudotime'],total.obs['latent_time']])
"""
array([[ 1.        , -0.0192098 ,  0.06366121,  0.06085888],
       [-0.0192098 ,  1.        ,  0.54726699,  0.54686801],
       [ 0.06366121,  0.54726699,  1.        ,  0.98841659],
       [ 0.06085888,  0.54686801,  0.98841659,  1.        ]])
"""

# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
( np.round(np.mean(v2s_mean),4) , np.round(np.mean(v2s_median),4) )

np.round(np.var(v2s_mean),4) # 
np.round(np.var(v2s_median),4) 
np.round(np.var(v2s_mean),4) # 0.0001
np.round(np.var(v2s_median),4) # 0.0299

# paired: 0.8597 0.9949
# shuffled: 0.0198 0.3261


####################################
# for algorithm evaluation purpose given expert annotated ground truth. 
# It contains a list of tuples in which stores the source cluster and target cluster of cells.
import unitvelo as utv
cluster_edges = [
    ("Blood progenitors 1", "Blood progenitors 2"),
    ("Blood progenitors 2", "Erythroid1"),
    ("Erythroid1", "Erythroid2"),
    ("Erythroid2", "Erythroid3")]
total_velo = total[:, total.var.loc[total.var['velocity_genes'] == True].index]
utv.evaluate(total_velo, cluster_edges, 'celltype', 'velocity')
# ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()