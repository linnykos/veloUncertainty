
import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

def utv_panINC_plots(split_seed):
    print(' ######################### seed='+str(split_seed)+': starts! ######################### ')
    dataset_long = 'pancreasINC'
    dataset_short = 'panINC'
    method = 'utv'
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

    total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
    # compute umap
    compute_umap(split1,dataset_short)
    compute_umap(split2,dataset_short)
    # add colors
    #raw = read_raw_adata(dataset_short)
    split1.uns['clusters_colors'] = total.uns['clusters_colors'].copy()
    split2.uns['clusters_colors'] = total.uns['clusters_colors'].copy()
    split1.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
    split2.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
    ## velocity
    print('seed='+str(split_seed)+': velocity plots')
    plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed)
    ## cosine similarity
    print('seed='+str(split_seed)+': cosine similarity')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    ## confidence
    print('seed='+str(split_seed)+': velocity confidence')
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)
    plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)
    print('correlation of velocity confidence between splits ='+str(np.round(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']),5)) )
    ## ptime
    print('seed='+str(split_seed)+': velocity pseudotime')
    scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)
    plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
    plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
    ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)
    ## latent time
    print('seed='+str(split_seed)+': latent time')
    scv.tl.recover_dynamics(split1,n_jobs=8)
    scv.tl.velocity_graph(split1,n_jobs=8)
    scv.tl.latent_time(split1)
    scv.tl.recover_dynamics(split2,n_jobs=8)
    scv.tl.velocity_graph(split2,n_jobs=8)
    scv.tl.latent_time(split2)
    plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)
    # shuffled cosine similarity
    print('seed='+str(split_seed)+': shuffled cosine similarity')
    v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
    print('mean of shuffled mean, mean of shuffled median, var of shuffled mean, var of shuffled median')
    print( np.round(np.mean(v2s_mean),4) , np.round(np.mean(v2s_median),4) )
    print(np.round(np.var(v2s_mean),4)) # 
    print(np.round(np.var(v2s_median),4)) 
    # intersected and unioned genes
    print('seed='+str(split_seed)+': cosine similarity using intersected genes')
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
    np.quantile(c1, [0.,.25,.5,.75,1.]) 
    print('seed='+str(split_seed)+': cosine similarity using unioned genes')
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    np.quantile(c2, [0.,.25,.5,.75,1.]) 
    print(' ######################### seed='+str(split_seed)+': all done! ######################### ')


utv_panINC_plots(323)
utv_panINC_plots(326)
utv_panINC_plots(329)



#total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
#split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
#split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
