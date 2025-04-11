import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import plot_velocity_scv_utv
from v4_functions import *

def plot_utv_ery_Mark(gene_set_name, split_seed, plot_total=True):
    dataset_long = 'erythroid'
    dataset_short = 'ery'
    method_prefix = 'utv'
    method = method_prefix + '_' + gene_set_name
    celltype_label = get_celltype_label(dataset_short)
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_outputAdded.h5ad')
    #total.obsm['X_umapOriginal'] = total.obsm['X_umap'].copy()
    total.obsm['X_pcaOriginal'] = total.obsm['X_pca'].copy()
    #del total.obsm['X_umap']
    del total.obsm['X_pca']
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded.h5ad')
    ## velocity
    if (plot_total):
        plot_velocity_scv_utv(adata_in=total,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
    plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
    plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
    print('####################### velocity plot done!')
    ## cosine similarity
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
    print('####################### cosine similarity plot done!')
    ## velocity confidence
    scv.tl.velocity_confidence(total)
    plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=celltype_label)
    print('####################### velocity confidence plot done!')

plot_total = True
for i in range(5):
    split_seed = [317, 320, 323, 326, 329][i]
    grid_seed = [227, 230, 233, 236, 239][i]
    gene_set_name = 'nMark' + str(grid_seed)
    plot_utv_ery_Mark(gene_set_name=gene_set_name, split_seed=split_seed, plot_total=plot_total)
    plot_utv_ery_Mark(gene_set_name=gene_set_name, split_seed=split_seed, plot_total=plot_total)
    plot_utv_ery_Mark(gene_set_name=gene_set_name, split_seed=split_seed, plot_total=plot_total)
    plot_total = False

for i in range(5):
    split_seed = [317, 320, 323, 326, 329][i]
    grid_seed = [227, 230, 233, 236, 239][i]
    gene_set_name = 'nMark' + str(grid_seed)
    plot_utv_ery_Mark(gene_set_name=gene_set_name, split_seed=split_seed, plot_total=True)
