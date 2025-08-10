dataset_long = 'greenleaf'
dataset_short = 'glf'

import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *


def create_plots(velo_method, split_seed, grid_seed):
    gene_set_name = 'nGPCgrid'+str(grid_seed)
    method = velo_method+'_'+gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'
    ## read data
    total = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_glf_'+method+'_total_v4.h5ad')
    split1 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
    split2 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
    total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
    ## compute UMAP
    compute_umap(total,dataset_short)
    compute_umap(split1,dataset_short)
    compute_umap(split2,dataset_short)
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
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

create_plots(velo_method='utv',split_seed=317, grid_seed=227)
print('################# plots of 317 done')

create_plots(velo_method='utv',split_seed=320, grid_seed=230)
print('################# plots of 320 done')

create_plots(velo_method='utv',split_seed=323, grid_seed=233)
print('################# plots of 323 done')

create_plots(velo_method='utv',split_seed=326, grid_seed=236)
print('################# plots of 326 done')

create_plots(velo_method='utv',split_seed=329, grid_seed=239)
print('################# plots of 329 done')

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'
def utv_create_cos_sim_df_across_seed(velo_method, split_seeds, grid_seeds):
    df = pd.DataFrame()
    for i in range(len(split_seeds)):
        split_seed = split_seeds[i]
        grid_seed = grid_seeds[i]
        gene_set_name = 'nGPCgrid'+str(grid_seed)
        method = velo_method+'_'+gene_set_name
        split1 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
        split2 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
        df['split'+str(split_seed)] = compute_cosine_similarity_union(split1,split2,method)[0]
        print(split_seed)
    return df

df = create_cos_sim_df_across_seed(velo_method='utv', split_seeds=[317,320,323,326,329], grid_seeds=[227,230,233,236,239])
df.to_csv(data_folder+'cos_sim_across_seeds_nGPCgrid_utv.csv')


print('################# all done')

