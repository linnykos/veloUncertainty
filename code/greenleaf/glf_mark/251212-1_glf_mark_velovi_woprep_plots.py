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

dataset_short = 'glf'
dataset_long = 'greenleaf'
velo_method = 'velovi_woprep'

######################################################
## compute umap

def compute_umap_GPC(adata):
    scv.pp.moments(adata, n_pcs=5, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=5) # used to be 10
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata) # since for nGPC genes, the neighbor structure might crash

# data_suffix = 'mark' or 'markctrl'
def glf_mark_velovi_create_plots(split_seed, method, data_suffix):
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/v4_greenleaf/seed'+str(split_seed)+'/'+method+'_'+data_suffix+'/'
    print('######## read data')
    total = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_total_'+dataset_short+'_'+data_suffix+'.h5ad')
    split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split1_'+dataset_short+'_'+data_suffix+'.h5ad') 
    split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split2_'+dataset_short+'_'+data_suffix+'.h5ad') 
    vae_total = VELOVI.load(data_folder+'vae_'+dataset_short+'_'+method+'_total_'+dataset_short+'_'+data_suffix+'.pt', total)
    vae_split1 = VELOVI.load(data_folder+'vae_'+dataset_short+'_'+method+'_split1_'+dataset_short+'_'+data_suffix+'.pt', split1)
    vae_split2 = VELOVI.load(data_folder+'vae_'+dataset_short+'_'+method+'_split2_'+dataset_short+'_'+data_suffix+'.pt', split2)
    total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
    ## add velovi outputs to adata
    print_message_with_time("############## Add velovi outputs to adata")
    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    add_velovi_outputs_to_adata(total, vae_total)
    compute_umap_GPC(split1)
    compute_umap_GPC(split2)
    compute_umap_GPC(total) # n_neighbors=30, n_pcs=40
    scv.tl.velocity_graph(total)
    scv.tl.velocity_graph(split1)
    scv.tl.velocity_graph(split2)
    ## save data
    total.write(data_folder+'adata_'+dataset_short+'_'+method+'_total_'+dataset_short+'_'+data_suffix+'_outputAdded.h5ad') # 
    print('######## write total to: '+data_folder+'adata_'+dataset_short+'_'+method+'_total_'+dataset_short+'_'+data_suffix+'_outputAdded.h5ad')
    split1.write(data_folder+'adata_'+dataset_short+'_'+method+'_split1_'+dataset_short+'_'+data_suffix+'_outputAdded.h5ad') # 
    print('######## write total to: '+data_folder+'adata_'+dataset_short+'_'+method+'_split1_'+dataset_short+'_'+data_suffix+'_outputAdded.h5ad')
    split2.write(data_folder+'adata_'+dataset_short+'_'+method+'_split2_'+dataset_short+'_'+data_suffix+'_outputAdded.h5ad') # 
    print('######## write total to: '+data_folder+'adata_'+dataset_short+'_'+method+'_split2_'+dataset_short+'_'+data_suffix+'_outputAdded.h5ad')
    ## plot velocity
    print('######## velocity')
    plot_velocity(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)
    ######################################################
    ## plot cosine similarity
    print('######## cosine similarity')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    ######################################################
    ## plot velo_conf
    print('######## velocity confidence')
    scv.tl.velocity_confidence(total)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)


#glf_mark_velovi_create_plots(split_seed=317, method=velo_method, data_suffix='mark')
#print('####################### 317 done')
glf_mark_velovi_create_plots(split_seed=320, method=velo_method, data_suffix='mark')
print('####################### 320 done')
glf_mark_velovi_create_plots(split_seed=323, method=velo_method, data_suffix='mark')
print('####################### 323 done')
glf_mark_velovi_create_plots(split_seed=326, method=velo_method, data_suffix='mark')
print('####################### 326 done')
glf_mark_velovi_create_plots(split_seed=329, method=velo_method, data_suffix='mark')
print('####################### 329 done')

print('####################### all done.')






