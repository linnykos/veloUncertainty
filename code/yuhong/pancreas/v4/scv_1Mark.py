import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd

from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'pancreas'
dataset_short = 'pan'
celltype_label = 'clusters'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
gene_set_prefix = 'Mark'

def run_scv_Mark(split_seed, grid_seed, gene_set_prefix):
    gene_set_name = gene_set_prefix+str(grid_seed)
    method = 'scv_'+gene_set_name
    print('################## Read data') 
    total = sc.read(data_folder+dataset_short+'_total_'+gene_set_name+'.h5ad')
    split1 = sc.read(data_folder+dataset_short+'_split1_'+gene_set_name+'.h5ad')
    split2 = sc.read(data_folder+dataset_short+'_split2_'+gene_set_name+'.h5ad')
    split1.layers['spliced_original'] = split1.layers['spliced'].copy()
    split1.layers['unspliced_original'] = split1.layers['unspliced'].copy()
    split2.layers['spliced_original'] = split1.layers['spliced'].copy()
    split2.layers['unspliced_original'] = split1.layers['unspliced'].copy()
    total.layers['spliced_original'] = total.layers['spliced'].copy()
    total.layers['unspliced_original'] = total.layers['unspliced'].copy()
    print('################## start velocity') 
    scv_compute_velocity(total, dataset_short) 
    scv_compute_velocity(split1, dataset_short) 
    scv_compute_velocity(split2, dataset_short) 
    print('################## end velocity') 
    total.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2.write_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')

## plots
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *


def plot_scv_Mark(split_seed, grid_seed, gene_set_prefix):
    gene_set_name = gene_set_prefix+str(grid_seed)
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


run_scv_Mark(split_seed=317, grid_seed=227, gene_set_prefix=gene_set_prefix)
plot_scv_Mark(split_seed=317, grid_seed=227, gene_set_prefix=gene_set_prefix)
print('##################### 317 done')

run_scv_Mark(split_seed=320, grid_seed=230, gene_set_prefix=gene_set_prefix)
plot_scv_Mark(split_seed=320, grid_seed=230, gene_set_prefix=gene_set_prefix)
print('##################### 320 done')

run_scv_Mark(split_seed=323, grid_seed=233, gene_set_prefix=gene_set_prefix)
plot_scv_Mark(split_seed=323, grid_seed=233, gene_set_prefix=gene_set_prefix)
print('##################### 323 done')

run_scv_Mark(split_seed=326, grid_seed=236, gene_set_prefix=gene_set_prefix)
plot_scv_Mark(split_seed=326, grid_seed=236, gene_set_prefix=gene_set_prefix)
print('##################### 326 done')

run_scv_Mark(split_seed=329, grid_seed=239, gene_set_prefix=gene_set_prefix)
plot_scv_Mark(split_seed=329, grid_seed=239, gene_set_prefix=gene_set_prefix)
print('##################### 329 done')


