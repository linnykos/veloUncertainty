split_seed = 320

import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_scv import *


def scv_larryMult_plots(split_seed):
    dataset_long = "larryMult"
    dataset_short = "larryMult"
    method = "scv"
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    adata_prefix = 'adata_'+dataset_short+'_'+method
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
    colors = ["#6e8ea1","#ffab6e","#dba8bc","#a0a0a0","#c4c88a","#87c3c9"]
    split1.uns['state_info_colors'] = colors
    split2.uns['state_info_colors'] = colors
    total.uns['state_info_colors'] = colors
    # gene correlation
    print('seed='+str(split_seed)+', gene correlation')
    plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)
    ## plot velocity
    print('seed='+str(split_seed)+', velocity of splits')
    plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)
    ## plot cosine similarity
    print('seed='+str(split_seed)+', cosine similarity plots')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    ## plot velocity confidence
    print('seed='+str(split_seed)+', velocity confidence plots')
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)
    print('seed='+str(split_seed)+', correlation of confidence btw splits='+str(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence'])))
    # pseudo time
    print('seed='+str(split_seed)+', pseudotime plots')
    scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)
    plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
    plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
    plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
    ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed,celltype_label='state_info')
    # latent time
    print('seed='+str(split_seed)+', latent time')
    scv.tl.latent_time(total)
    scv.tl.latent_time(split1)
    scv.tl.latent_time(split2)
    plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    # shuffled cosine similarity
    print('seed='+str(split_seed)+', shuffled cosine similarity')
    v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
    print('shuffled mean of mean='+str(np.round(np.mean(v2s_mean),4) ))
    print('shuffled mean of median='+str(np.round(np.mean(v2s_median),4) ))
    print('shuffled variance of mean'+str(np.round(np.var(v2s_mean),4))) 
    print('shuffled variance of median'+str(np.round(np.var(v2s_median),4))) 
    # intersection vs union
    print('seed='+str(split_seed)+', cosine similarity using intersected vs unioned genes')
    print('cosine similarity of intersected genes')
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) 
    print(np.quantile(c1,[0.,.25,.5,.75,1.])) 
    print('cosine similarity of unioned genes')
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    print(np.quantile(c2,[0.,.25,.5,.75,1.]))
    print('All done for split_seed='+str(split_seed))



scv_larryMult_plots(split_seed=320)
scv_larryMult_plots(split_seed=323)
scv_larryMult_plots(split_seed=326)
scv_larryMult_plots(split_seed=329)

