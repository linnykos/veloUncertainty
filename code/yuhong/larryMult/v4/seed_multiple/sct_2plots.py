import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import torch
import random
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *
from v4_functions import *


def sct_larryMult_plots(split_seed,timesteps=[i/50 for i in range(1,11)]):
    method = 'sct'
    dataset_long = 'larryMult'
    dataset_short = 'larryMult'
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
    #timestep = 0.61
    #df = pd.read_csv(data_folder+dataset_short+'_sct_velo_timestep.csv')
    #timestep = df.iloc[np.where(df['pseudotime_corr']==np.max(df['pseudotime_corr']))[0][0]]['time'] # 0.61
    adata_prefix = 'adata_'+dataset_short+'_'+method
    tnode_prefix = 'tnode_'+dataset_short+'_'+method
    total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
    tnode_total = sct.predict.load_model('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed317/sct/'+tnode_prefix+'_total_v4.pth')
    tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
    tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')
    colors = ["#6e8ea1","#ffab6e","#dba8bc","#a0a0a0","#c4c88a","#87c3c9"]
    split1.uns['state_info_colors'] = colors
    split2.uns['state_info_colors'] = colors
    total.uns['state_info_colors'] = colors
    compute_umap(split1, dataset_short)
    compute_umap(split2, dataset_short)
    compute_umap(total, dataset_short)
    total.obsm['X_umapOriginal'] = total.obsm['X_umap'].copy()
    total.obsm['X_umapOriginal'][:,0] = np.array(total.obs['SPRING-x'])
    total.obsm['X_umapOriginal'][:,1] = np.array(total.obs['SPRING-y'])
    split1.obsm['X_umapOriginal'] = split1.obsm['X_umap'].copy()
    split1.obsm['X_umapOriginal'][:,0] = np.array(split1.obs['SPRING-x'])
    split1.obsm['X_umapOriginal'][:,1] = np.array(split1.obs['SPRING-y'])
    split2.obsm['X_umapOriginal'] = split2.obsm['X_umap'].copy()
    split2.obsm['X_umapOriginal'][:,0] = np.array(split2.obs['SPRING-x'])
    split2.obsm['X_umapOriginal'][:,1] = np.array(split2.obs['SPRING-y'])
    sct_seed=615
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    total.layers['velocity'] = compute_sct_avg_velocity(tnode_total, timesteps)
    split1.layers['velocity'] = compute_sct_avg_velocity(tnode_split1, timesteps) 
    split2.layers['velocity'] = compute_sct_avg_velocity(tnode_split2, timesteps) 
    print('seed'+str(split_seed)+': write outputAdded adata')
    split1.write_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
    split2.write_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
    print('seed'+str(split_seed)+': plot gene correlations')
    plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)
    # vector field
    print('seed'+str(split_seed)+': plot vector field')
    plot_vf_umap(adata_in=split1, data_version="split1",data=dataset_short,method=method,fig_folder=fig_folder)
    plot_vf_umap(adata_in=split2, data_version="split2",data=dataset_short,method=method,fig_folder=fig_folder)
    ptime_sct_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)
    ## velocity
    print('seed'+str(split_seed)+': plot velocity')
    plot_sct_velocity(adata_in=split1,data_version='split1',dataset=dataset_short,fig_folder=fig_folder,recompute=True,method='sct',celltype_label=None)
    plot_sct_velocity(adata_in=split2,data_version='split2',dataset=dataset_short,fig_folder=fig_folder,recompute=True,method='sct',celltype_label=None)
    ## cosine similarity
    print('seed'+str(split_seed)+': plot cosine similarity')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    ## plot velo_conf
    print('seed'+str(split_seed)+': plot velocity confidence')
    if (not 'velocity_confidence' in total.obs.columns):
        scv.tl.velocity_confidence(total)
        scv.tl.velocity_confidence(split1)
        scv.tl.velocity_confidence(split2)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)
    print('correlation of velocity confidence between splits='+str(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence'])))
    ## ptime
    print('seed'+str(split_seed)+': plot pseudotime')
    scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)
    plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
    plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
    plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
    ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)
    # latent time
    print('seed'+str(split_seed)+': plot latent time')
    if not 'latent_time' in split1.obs.columns:
        scv.tl.recover_dynamics(total,n_jobs=8)
        scv.tl.latent_time(total)
        plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
        scv.tl.recover_dynamics(split1,n_jobs=8)
        scv.tl.latent_time(split1)
        plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
        scv.tl.recover_dynamics(split2,n_jobs=8)
        scv.tl.latent_time(split2)
        plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    print('seed'+str(split_seed)+': shuffled cosine similarity')
    v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
    print('mean of shuffled means'+str(np.round(np.mean(v2s_mean),5)))
    print('mean of shuffled medians'+str(np.round(np.mean(v2s_median),5)))
    print('variance of shuffled mean'+str(np.round(np.var(v2s_mean),4)))
    print('variance of shuffled median'+str(np.round(np.var(v2s_median),4)))
    print('seed='+str(split_seed)+': Compute cosine similarity using intersected genes and unioned genes')
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
    print(np.quantile(c1,[0.,.25,.5,.75,1.]) )
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    print(np.quantile(c2,[0.,.25,.5,.75,1.]) )
    print('############## seed='+str(split_seed)+': All done')


sct_larryMult_plots(split_seed=320)
sct_larryMult_plots(split_seed=323)
sct_larryMult_plots(split_seed=326)
sct_larryMult_plots(split_seed=329)

exit()

def compute_sct_avg_velocity(tnode,timesteps):
    v_shape = tnode.adata.shape
    v = np.zeros(v_shape)
    for t in timesteps:
        v += compute_sctour_velocity(tnode, timestep=t)
    return v/len(timesteps)

def sct_larryMult_write_splits_outputAdded(split_seed):
    method = 'sct'
    dataset_long = 'larryMult'
    dataset_short = 'larryMult'
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
    #adata_prefix = 'adata_'+dataset_short+'_'+method
    tnode_prefix = 'tnode_'+dataset_short+'_'+method
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
    tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
    tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')
    colors = ["#6e8ea1","#ffab6e","#dba8bc","#a0a0a0","#c4c88a","#87c3c9"]
    split1.uns['state_info_colors'] = colors
    split2.uns['state_info_colors'] = colors
    compute_umap(split1, dataset_short)
    compute_umap(split2, dataset_short)
    split1.obsm['X_umapOriginal'] = split1.obsm['X_umap'].copy()
    split1.obsm['X_umapOriginal'][:,0] = np.array(split1.obs['SPRING-x'])
    split1.obsm['X_umapOriginal'][:,1] = np.array(split1.obs['SPRING-y'])
    split2.obsm['X_umapOriginal'] = split2.obsm['X_umap'].copy()
    split2.obsm['X_umapOriginal'][:,0] = np.array(split2.obs['SPRING-x'])
    split2.obsm['X_umapOriginal'][:,1] = np.array(split2.obs['SPRING-y'])
    # compute avg velocity
    sct_seed=615
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    timesteps=[i/50 for i in range(1,11)]
    #total.layers['velocity'] = compute_sct_avg_velocity(tnode_total, timesteps)
    split1.layers['velocity'] = compute_sct_avg_velocity(tnode_split1, timesteps) 
    split2.layers['velocity'] = compute_sct_avg_velocity(tnode_split2, timesteps)
    #split1.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=timestep) 
    #split2.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=timestep) 
    split1.write_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
    split2.write_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')

#sct_larryMult_write_splits_outputAdded(split_seed=320)
#sct_larryMult_write_splits_outputAdded(split_seed=323)
#sct_larryMult_write_splits_outputAdded(split_seed=326)
#sct_larryMult_write_splits_outputAdded(split_seed=329)


"""
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
tnode_total = sct.predict.load_model(data_folder+tnode_prefix+'_total_v4.pth')
tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')

colors = ["#6e8ea1","#ffab6e","#dba8bc","#a0a0a0","#c4c88a","#87c3c9"]
split1.uns['state_info_colors'] = colors
split2.uns['state_info_colors'] = colors
total.uns['state_info_colors'] = colors

compute_umap(split1, dataset_short)
compute_umap(split2, dataset_short)
compute_umap(total, dataset_short)

total.obsm['X_umapOriginal'] = total.obsm['X_umap'].copy()
total.obsm['X_umapOriginal'][:,0] = np.array(total.obs['SPRING-x'])
total.obsm['X_umapOriginal'][:,1] = np.array(total.obs['SPRING-y'])

split1.obsm['X_umapOriginal'] = split1.obsm['X_umap'].copy()
split1.obsm['X_umapOriginal'][:,0] = np.array(split1.obs['SPRING-x'])
split1.obsm['X_umapOriginal'][:,1] = np.array(split1.obs['SPRING-y'])

split2.obsm['X_umapOriginal'] = split2.obsm['X_umap'].copy()
split2.obsm['X_umapOriginal'][:,0] = np.array(split2.obs['SPRING-x'])
split2.obsm['X_umapOriginal'][:,1] = np.array(split2.obs['SPRING-y'])

total.layers['velocity'] = compute_sctour_velocity(tnode_total, timestep=timestep) #+
split1.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=timestep) #-
split2.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=timestep) #+

total.write(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1.write(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2.write(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)

"""



