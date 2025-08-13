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

def compute_latent_time(adata):
    scv.tl.recover_dynamics(adata,n_jobs=8)
    scv.tl.velocity_graph(adata,n_jobs=8)
    scv.tl.latent_time(adata)

def velovi_panINC_plots(split_seed):
    print_message_with_time('####################### seed'+str(split_seed)+': starts! #######################') 
    method = 'velovi'
    dataset_short = 'panINC'
    dataset_long = 'pancreasINC'
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    print_message_with_time(str(split_seed)+': read data') 
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
    vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split1_v4.pt', split1)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
    vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split2_v4.pt', split2)
    total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
    ## add velovi outputs to adata
    print_message_with_time(str(split_seed)+': add velovi output to splits') 
    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    ## compute umap
    print_message_with_time(str(split_seed)+': compute umap') 
    print_message_with_time("Compute umap")
    compute_umap(split1, dataset_short)
    compute_umap(split2, dataset_short)
    ## write data
    split1.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
    split2.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
    # gene correlation
    print_message_with_time(str(split_seed)+': gene correlation') 
    plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)
    ## plot velocity
    print_message_with_time(str(split_seed)+': velocity plots') 
    plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)
    ## plot cosine similarity
    print_message_with_time(str(split_seed)+': cosine similarity plots') 
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    ## plot velo_conf
    print_message_with_time(str(split_seed)+': velocity confidence plots') 
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)
    print('correlation of velocity confidence between splits='+str(np.round(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']),5)))
    # shuffled cosine similarity
    print_message_with_time(str(split_seed)+': shuffled cosine similarity')
    v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
    print('mean of shuffled mean, mean of shuffled median, variance of shuffled mean, variance of shuffled median')
    print(np.round(np.mean(v2s_mean),4))
    print(np.round(np.mean(v2s_median),4))
    print(np.round(np.var(v2s_mean),4)) 
    print(np.round(np.var(v2s_median),4)) 
    # cosine similarity using intersected and unioned genes
    print_message_with_time(str(split_seed)+': cosine similarity of intersected genes')
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
    print(np.quantile(c1,[0.,.25,.5,.75,1.]) )
    print_message_with_time(str(split_seed)+': cosine similarity of unioned genes')
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    print(np.quantile(c2,[0.,.25,.5,.75,1.]) )
    print_message_with_time('####################### seed'+str(split_seed)+': all done! #######################') 

velovi_panINC_plots(split_seed=320)
velovi_panINC_plots(split_seed=323)
velovi_panINC_plots(split_seed=326)
velovi_panINC_plots(split_seed=329)
