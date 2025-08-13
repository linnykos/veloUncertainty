dataset_long = 'pancreas'
dataset_short = 'pan'
method = 'scv'
split_seed=317

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

celltype_label = get_celltype_label(dataset_short)

total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_total_v4.h5ad')
split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split1_v4.h5ad')
split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split2_v4.h5ad')

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
# 716 1102
np.quantile(c1,[0.,.25,.5,.75,1.]) # [-0.96536275,  0.18747088,  0.54379243,  0.8632963 ,  0.99629776]
np.quantile(c2,[0.,.25,.5,.75,1.]) # [-0.96297644,  0.16518448,  0.4658463 ,  0.83853388,  0.99470953]

## confidence
scv.tl.velocity_confidence(total)
scv.tl.velocity_confidence(split1)
scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)

###################################
# method-selected gene corr
plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)

####################################
# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
np.round(np.mean(v2s_mean),5) # 
np.round(np.mean(v2s_median),5) # 

# paired: 0.492 0.4658
# shuffled: 0.41013 0.2769
