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

def velovi_larryMult_plots(split_seed):
    method = 'velovi'
    dataset_short = 'larryMult'
    dataset_long = 'larryMult'
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
    vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split1_v4.pt', split1)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
    vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split2_v4.pt', split2)
    total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
    vae_total = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed317/'+method+'/vae_'+dataset_short+'_'+method+'_total_v4.pt', total)
    colors = ['#6e8ea1','#ffab6e','#dba8bc','#a0a0a0','#c4c88a','#87c3c9']
    split1.uns['state_info_colors'] = colors
    split2.uns['state_info_colors'] = colors
    total.uns['state_info_colors'] = colors
    ## add velovi outputs to adata
    print_message_with_time('############## Add velovi outputs to adata')
    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    add_velovi_outputs_to_adata(total, vae_total)
    ## compute umap
    print_message_with_time('############## Compute umap')
    compute_umap(split1, dataset_short)
    compute_umap(split2, dataset_short)
    compute_umap(total, dataset_short)
    ## add umapOriginal
    print_message_with_time('############## Add umapOriginal to adata')
    total.obsm['X_umapOriginal'] = total.obsm['X_umap'].copy()
    total.obsm['X_umapOriginal'][:,0] = np.array(total.obs['SPRING-x'])
    total.obsm['X_umapOriginal'][:,1] = np.array(total.obs['SPRING-y'])
    split1.obsm['X_umapOriginal'] = split1.obsm['X_umap'].copy()
    split1.obsm['X_umapOriginal'][:,0] = np.array(split1.obs['SPRING-x'])
    split1.obsm['X_umapOriginal'][:,1] = np.array(split1.obs['SPRING-y'])
    split2.obsm['X_umapOriginal'] = split2.obsm['X_umap'].copy()
    split2.obsm['X_umapOriginal'][:,0] = np.array(split2.obs['SPRING-x'])
    split2.obsm['X_umapOriginal'][:,1] = np.array(split2.obs['SPRING-y'])
    # gene correlation
    print_message_with_time('############## seed='+str(split_seed)+': Plot gene correlation')
    plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)
    ## plot velocity
    print_message_with_time('############## seed='+str(split_seed)+': Add velovi outputs to adata')
    plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed)
    ## plot cosine similarity
    print_message_with_time('############## seed='+str(split_seed)+': Plot cosine similarity')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    print_message_with_time('############## seed='+str(split_seed)+': Plot velocity confidence')
    if (not 'velocity_confidence' in total.obs.columns):
        scv.tl.velocity_confidence(total)
        scv.tl.velocity_confidence(split1)
        scv.tl.velocity_confidence(split2)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)
    print('correlation of velocity confidence between splits='+str(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence'])))
    ## ptime
    print_message_with_time('############## seed='+str(split_seed)+': Plot pseudotime')
    if not 'velocity_pseudotime' in split1.obs.columns:
        scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
        scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
        scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)
    plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
    plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
    ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed,celltype_label='state_info')
    # shuffled cosine similarity
    print_message_with_time('############## seed='+str(split_seed)+': Compute shuffled cosine similarity')
    v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
    print('mean of shuffled mean'+str(np.round(np.mean(v2s_mean),4)))
    print('mean of shuffled median'+str(np.round(np.mean(v2s_median),4)))
    print('variance of shuffled mean'+str(np.round(np.var(v2s_mean),4)))
    print('variance of shuffled median'+str(np.round(np.var(v2s_median),4)))
    print_message_with_time('############## seed='+str(split_seed)+': Compute cosine similarity using intersected genes and unioned genes')
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
    print(np.quantile(c1,[0.,.25,.5,.75,1.]) )
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    print(np.quantile(c2,[0.,.25,.5,.75,1.]) )
    print_message_with_time('############## seed='+str(split_seed)+': All done')


velovi_larryMult_plots(split_seed=320)
velovi_larryMult_plots(split_seed=323)
velovi_larryMult_plots(split_seed=326)
velovi_larryMult_plots(split_seed=329)

######################################################

'''
if not 'latent_time' in split1.obs.columns:
    scv.tl.latent_time(total)
    scv.tl.latent_time(split1)
    scv.tl.latent_time(split2)

# KeyError: 'fit_likelihood'
plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
'''


# paired: 0.1801 0.1876
# shuffled: 0.1024 0.0901


'''
## write data
split1.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
split2.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
total.write_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad')

split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)

vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/'+method+'/vae_'+dataset_short+'_'+method+'_split1_v4.pt', split1)
vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/'+method+'/vae_'+dataset_short+'_'+method+'_split2_v4.pt', split2)
vae_total = VELOVI.load(data_folder+'v4_'+dataset_long+'/'+method+'/vae_'+dataset_short+'_'+method+'_total_v4.pt', total)
'''