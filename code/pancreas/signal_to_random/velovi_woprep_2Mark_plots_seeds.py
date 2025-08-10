dataset_short = 'pan'
dataset_long = 'pancreas'
gene_set_prefix = 'Mark'
method_prefix = 'velovi_woprep'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_velovi import *

######################################################
## compute umap
def compute_umap_GPC(adata):
    scv.pp.moments(adata, n_pcs=5, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=5) # used to be 10
    sc.tl.umap(adata)
    #scv.tl.velocity_graph(adata) # since for nGPC genes, the neighbor structure might crash


def velovi_create_plots(split_seed, grid_seed, velo_method):
    gene_set_name = gene_set_prefix+str(grid_seed)
    method = velo_method+'_'+gene_set_name
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    print('######## read data')
    total = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_total_GPC.h5ad') # 
    split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad') # 
    split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad') # 
    vae_total = VELOVI.load(data_folder+'vae_'+dataset_short+'_'+method+'_total_GPC.pt', total)
    vae_split1 = VELOVI.load(data_folder+'vae_'+dataset_short+'_'+method+'_split1_GPC.pt', split1)
    vae_split2 = VELOVI.load(data_folder+'vae_'+dataset_short+'_'+method+'_split2_GPC.pt', split2)
    ## add velovi outputs to adata
    print_message_with_time("############## Add velovi outputs to adata")
    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    add_velovi_outputs_to_adata(total, vae_total)
    print(method)
    compute_umap_GPC(split1)
    compute_umap_GPC(split2)
    compute_umap_GPC(total) # n_neighbors=30, n_pcs=40
    scv.tl.velocity_graph(total)
    total.write(data_folder+'adata_'+dataset_short+'_'+method+'_total_GPC_outputAdded.h5ad') # 
    split1.write(data_folder+'adata_'+dataset_short+'_'+method+'_split1_GPC_outputAdded.h5ad') # 
    split2.write(data_folder+'adata_'+dataset_short+'_'+method+'_split2_GPC_outputAdded.h5ad') # 
    ## plot velocity
    print('######## velocity')
    plot_velocity(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
    #plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
    #plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)
    ######################################################
    ## plot cosine similarity
    print('######## cosine similarity')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) 
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    print( np.quantile(c1,[0.,.25,.5,.75,1.]) )
    print( np.quantile(c2,[0.,.25,.5,.75,1.]) )
    ######################################################
    ## plot velo_conf
    print('######## velocity confidence')
    if (not 'velocity_confidence' in total.obs.columns):
        scv.tl.velocity_confidence(total)
    #plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)


velovi_create_plots(split_seed=317, grid_seed=227, velo_method=method_prefix)
print('####################### 317 done')
velovi_create_plots(split_seed=320, grid_seed=230, velo_method=method_prefix)
print('####################### 320 done')
velovi_create_plots(split_seed=323, grid_seed=233, velo_method=method_prefix)
print('####################### 323 done')
velovi_create_plots(split_seed=326, grid_seed=236, velo_method=method_prefix)
print('####################### 326 done')
velovi_create_plots(split_seed=329, grid_seed=239, velo_method=method_prefix)
print('####################### 329 done')

print('####################### all done.')