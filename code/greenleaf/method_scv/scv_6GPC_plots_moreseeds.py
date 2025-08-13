dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'scv_GPC'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

celltype_label = 'cluster_name'
total = sc.read(data_folder+'seed317/scv_GPC/adata_'+dataset_short+'_'+method+'_total_GPC.h5ad')
total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()

for split_seed in [320, 323, 326, 329]:
    print('######################## seed'+str(split_seed)+' starts!')
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    split1 = sc.read(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
    split2 = sc.read(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')
    split1.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
    split2.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
    ## velocity
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
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)
    plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=celltype_label)

    ## ptime
    print('############### pseudo-time')
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)
    plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
    plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
    ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)

    print('############### latent-time')
    scv.tl.latent_time(split1)
    plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    scv.tl.latent_time(split2)
    plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)
    print('################ seed'+str(split_seed)+' done!')