import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'scv_nGPC227'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
#fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

print('################## Read data') 
total = sc.read(data_folder+'glf_total_nGPC227.h5ad')
split1 = sc.read(data_folder+'seed317_glf_split1_nGPC227.h5ad')
split2 = sc.read(data_folder+'seed317_glf_split2_nGPC227.h5ad')

split1.layers['spliced_original'] = split1.layers['spliced'].copy()
split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()
split2.layers['spliced_original'] = split1.layers['spliced'].copy()
split2.layers['unspliced_original'] = split1.layers['unspliced'].copy()
total.layers['spliced_original'] = total.layers['spliced'].copy()
total.layers['unspliced_original'] = total.layers['unspliced'].copy()

total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
split1.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
split2.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()


print('################## start velocity') 
scv_compute_velocity(total,dataset_short) 
scv_compute_velocity(split1,dataset_short) 
scv_compute_velocity(split2,dataset_short) 
print('################## end velocity') 

total.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_GPC.h5ad')
split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')

## plots
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

celltype_label = 'cluster_name'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_GPC.h5ad')
split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')

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
print('intersected and unioned gene cosine similarity quantiles:\n')
np.quantile(c1,[0.,.25,.5,.75,1.]) 
np.quantile(c2,[0.,.25,.5,.75,1.]) 
print(np.corrcoef(c1,c2))

## confidence
print('############### velocity confidence')
scv.tl.velocity_confidence(total)
scv.tl.velocity_confidence(split1)
scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=celltype_label)
