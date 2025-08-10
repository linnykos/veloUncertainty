import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

def run_scv_nGPCblk(split_seed, grid_seed):
    gene_set_name = 'nGPCblk'+str(grid_seed)
    method = 'scv_'+gene_set_name
    print('################## Read data') 
    total = sc.read(data_folder+'glf_total_'+gene_set_name+'.h5ad')
    split1 = sc.read(data_folder+'glf_seed'+str(split_seed)+'_split1_'+gene_set_name+'.h5ad')
    split2 = sc.read(data_folder+'glf_seed'+str(split_seed)+'_split2_'+gene_set_name+'.h5ad')

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
    scv_compute_velocity(total, dataset_short) 
    scv_compute_velocity(split1, dataset_short) 
    scv_compute_velocity(split2, dataset_short) 
    print('################## end velocity') 

    total.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')

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

def plot_scv_nGPCblk(split_seed, grid_seed):
    gene_set_name = 'nGPCblk'+str(grid_seed)
    method = 'scv_'+gene_set_name
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

    total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
    print('############### read data')
    ## velocity
    plot_velocity_scv_utv(adata_in=total,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,
                          split_seed=split_seed,celltype_label=celltype_label)
    plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,
                          split_seed=split_seed,celltype_label=celltype_label)
    plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,
                          split_seed=split_seed,celltype_label=celltype_label)
    print('############### velocity plotted')
    ## cosine similarity
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,
                           fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,
                                            method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,
                                               method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    print('intersected and unioned gene cosine similarity quantiles:\n')
    print( np.quantile(c1,[0.,.25,.5,.75,1.]) )
    print( np.quantile(c2,[0.,.25,.5,.75,1.]) )
    print('############### replicate coherence done')
    ## confidence
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

    plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=celltype_label)
    print('############### velocity confidence done')


run_scv_nGPCblk(split_seed=317, grid_seed=227)
plot_scv_nGPCblk(split_seed=317, grid_seed=227)
print('#################### seed317 done')

run_scv_nGPCblk(split_seed=320, grid_seed=230)
plot_scv_nGPCblk(split_seed=320, grid_seed=230)
print('#################### seed320 done')

run_scv_nGPCblk(split_seed=323, grid_seed=233)
plot_scv_nGPCblk(split_seed=323, grid_seed=233)
print('#################### seed323 done')

run_scv_nGPCblk(split_seed=326, grid_seed=236)
plot_scv_nGPCblk(split_seed=326, grid_seed=236)
print('#################### seed326 done')

run_scv_nGPCblk(split_seed=329, grid_seed=239)
plot_scv_nGPCblk(split_seed=329, grid_seed=239)
print('#################### seed329 done')


def scv_create_cos_sim_df_across_seed(split_seeds, grid_seeds):
    velo_method = 'scv'
    df = pd.DataFrame()
    for i in range(len(split_seeds)):
        split_seed = split_seeds[i]
        grid_seed = grid_seeds[i]
        gene_set_name = 'nGPCblk'+str(grid_seed)
        method = velo_method+'_'+gene_set_name
        split1 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
        split2 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
        df['split'+str(split_seed)] = compute_cosine_similarity_union(split1,split2,method)[0]
        print(split_seed)
    return df

df = scv_create_cos_sim_df_across_seed(split_seeds=[317,320,323,326,329], grid_seeds=[227,230,233,236,239])
df.to_csv(data_folder+'cos_sim_across_seeds_nGPCblk_scv.csv')
