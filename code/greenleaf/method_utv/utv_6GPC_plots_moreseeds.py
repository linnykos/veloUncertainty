dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'utv_GPC'

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

total = sc.read(data_folder+'seed317/utv/adata_glf_utv_total_v4_outputAdded.h5ad')
total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
compute_umap(total,dataset_short)

for split_seed in [320,323,326,329]:
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    split1 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/seed'+str(split_seed)+'_glf_split1_GPC_utv.h5ad')
    split2 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/seed'+str(split_seed)+'_glf_split2_GPC_utv.h5ad')
    compute_umap(split1,dataset_short)
    compute_umap(split2,dataset_short)

    ## velocity
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
    print(np.quantile(c1, [0.,.25,.5,.75,1.]) )
    print(np.quantile(c2, [0.,.25,.5,.75,1.]) )

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

