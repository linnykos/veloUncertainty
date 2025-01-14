dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'utv_GPC'
split_seed=317

import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total = sc.read(data_folder+'seed'+str(split_seed)+'/utv/glf_total_GPC_utv.h5ad')
split1 = sc.read(data_folder+'seed'+str(split_seed)+'/utv/seed'+str(split_seed)+'_glf_split1_GPC_utv.h5ad')
split2 = sc.read(data_folder+'seed'+str(split_seed)+'/utv/seed'+str(split_seed)+'_glf_split2_GPC_utv.h5ad')
total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()

compute_umap(total,dataset_short)
compute_umap(split1,dataset_short)
compute_umap(split2,dataset_short)

split1.obs['split'] = 'split1'
split2.obs['split'] = 'split2'

sc.tl.pca(split1)
sc.tl.pca(split2)
import matplotlib.pyplot as plt
plt.figure(figsize=(4, 4))
plt.scatter(split1.obsm['X_pca'][:, 0], split1.obsm['X_pca'][:, 1], s=10, label='split1', alpha=0.3, color='blue')
plt.scatter(split2.obsm['X_pca'][:, 0], split2.obsm['X_pca'][:, 1], s=10, label='split2', alpha=0.3, color='orange')
plt.legend()
plt.title('PCA of split1 and 2')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.savefig(fig_folder+'utv_GPC_pca_overlay.png', dpi=100)

plt.figure(figsize=(4, 4))
plt.scatter(split1.obsm['X_umap'][:, 0], split1.obsm['X_umap'][:, 1], s=10, label='split1', alpha=0.3, color='blue')
plt.scatter(split2.obsm['X_umap'][:, 0], split2.obsm['X_umap'][:, 1], s=10, label='split2', alpha=0.3, color='orange')
plt.legend()
plt.title('UMAP of split1 and 2')
plt.savefig(fig_folder+'utv_GPC_umap_overlay.png', dpi=100)

## velocity
plot_velocity_scv_utv(adata_in=total,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed)

## cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

print('cosine similarity using gene intersection and union')
c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)

np.quantile(c1, [0.,.25,.5,.75,1.]) 
np.quantile(c2, [0.,.25,.5,.75,1.]) 


## confidence
if (not 'velocity_confidence' in split1.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)


print('shuffled cosine similarity')
# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
print('shuffled mean of mean and median: ')
( np.round(np.mean(v2s_mean),4) , np.round(np.mean(v2s_median),4) )

