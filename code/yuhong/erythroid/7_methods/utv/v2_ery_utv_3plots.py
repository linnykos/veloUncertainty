import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import bbknn

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'utv'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")

######################################################
## compute umap
def utv_compute_umap(adata):
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("************ batch correction done ************")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # used to be 10
    sc.tl.umap(adata)

utv_compute_umap(total)
utv_compute_umap(split1)
utv_compute_umap(split2)

######################################################
## plot velocity - done in 1data.py, copied (and modified variable names here)
plot_velocity_scv_utv(adata_in=total,adata_raw=raw,fig_folder=fig_folder,fig_info='total',dataset=dataset_short,method=method)
plot_velocity_scv_utv(adata_in=split1,adata_raw=raw,fig_folder=fig_folder,fig_info='split1',dataset=dataset_short,method=method)
plot_velocity_scv_utv(adata_in=split2,adata_raw=raw,fig_folder=fig_folder,fig_info='split2',dataset=dataset_short,method=method)

plot_velocity_scv_utv(adata_in=total,adata_raw=raw,fig_folder=fig_folder,fig_info='recompF_total',dataset=dataset_short,method=method,recompute=False)
plot_velocity_scv_utv(adata_in=split1,adata_raw=raw,fig_folder=fig_folder,fig_info='recompF_split1',dataset=dataset_short,method=method,recompute=False)
plot_velocity_scv_utv(adata_in=split2,adata_raw=raw,fig_folder=fig_folder,fig_info='recompF_split2',dataset=dataset_short,method=method,recompute=False)

######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)

plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_hist_by_celltype(split1,split2,total,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## plot velo_conf
plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_veloConf_hist(adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,text_x=.7)

######################################################
## ptime
plot_pseudotime(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,method=method,fig_folder=fig_folder)

ptime_correlation_scatter_plot(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder)

# Spearman's corr
ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder,time_label='velocity_pseudotime')
ptime_correlation_scatter_spearman(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder,time_label='velocity_pseudotime')
ptime_correlation_scatter_spearman(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder,time_label='velocity_pseudotime')

# latent time
if not 'latent_time' in split1.obs.columns:
    scv.tl.latent_time(total)
    scv.tl.latent_time(split1)
    scv.tl.latent_time(split2)

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder,time_label='latent_time')
ptime_correlation_scatter_spearman(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder,time_label='latent_time')
ptime_correlation_scatter_spearman(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder,time_label='latent_time')


exit()
######################################################
## plot velocity - done in 1data.py, copied (and modified variable names here)
raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
scv.pl.velocity_embedding_stream(total, basis='umap',color="sequencing.batch",save=fig_folder+"velocity/"+dataset_short+"_"+method+"_total_umapCompute_byBatch.png")

def plot_velocities_scv_utv(adata_in,adata_raw,fig_info,dataset,method,color_label=None):
    if dataset=="ery":
        color_label = 'celltype'
    elif dataset=="pan":
        color_label = 'clusters'
    data_method = dataset+"_"+method
    # umapCompute
    adata=adata_in.copy()
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=color_label,save=fig_folder+"velocity/"+data_method+"_"+fig_info+"_umapCompute.png")
    # umapOriginal
    adata=adata_in.copy()
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=color_label,save=fig_folder+"velocity/"+data_method+"_"+fig_info+"_umapOriginal.png")    
        
plot_velocities_scv_utv(adata_in=total,adata_raw=raw,fig_info="total",dataset=dataset_short,method=method)
plot_velocities_scv_utv(adata_in=split1,adata_raw=raw,fig_info="split1",dataset=dataset_short,method=method)
plot_velocities_scv_utv(adata_in=split2,adata_raw=raw,fig_info="split2",dataset=dataset_short,method=method)

######################################################
## plot cosine similarity
velo_genes_split1 = split1.var.index#[~np.isnan(split1.layers['velocity'][0])]
velo_genes_split2 = split2.var.index#[~np.isnan(split2.layers['velocity'][0])]
common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0])) # 1481

velo_df1 = pd.DataFrame(split1.layers['velocity'], columns=split1.var.index.tolist())
velo_df2 = pd.DataFrame(split2.layers['velocity'], columns=split2.var.index.tolist())
cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))
np.quantile(cos_sim,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

# same as compute_cosine_similarity_scv actually
def compute_cosine_similarity_utv(adata_split1,adata_split2):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0])) 
    velo_df1 = pd.DataFrame(adata_split1.layers['velocity'], columns=adata_split1.var.index.tolist())
    velo_df2 = pd.DataFrame(adata_split2.layers['velocity'], columns=adata_split2.var.index.tolist())
    cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))
    return cos_sim, common_genes_velocity.shape[0] # return cosine similarity and number of common genes in velocity computation

compute_cosine_similarity_utv(split1,split2)

def compute_cosine_similarity(adata_split1,adata_split2,method):
    if method=="scv":
        return compute_cosine_similarity_scv(adata_split1,adata_split2)
    elif method=="utv":
        return compute_cosine_similarity_utv(adata_split1,adata_split2)

import matplotlib.pyplot as plt
# cosine similarity
def plot_cosine_similarity(adata_split1,adata_split2,adata_total,adata_raw,dataset,method,fig_folder=fig_folder,text_x=None,text_y=None):
    cos_sim, Ngenes = compute_cosine_similarity(adata_split1,adata_split2,method)
    dataset_method = dataset+'_'+method
    # histogram
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='dimgray',color='powderblue') 
    max_frequency = np.max(counts)
    if text_x is None:
        text_x = np.quantile(cos_sim,[.05])[0]
    if text_y is None:
        text_y = max_frequency/2
    plt.axvline(np.mean(cos_sim), color='salmon', linestyle='dashed', linewidth=1.5) ## add mean
    plt.text(text_x,text_y,'mean='+str(np.round(np.mean(cos_sim),4))+', median='+str(np.round(np.median(cos_sim),4)), color='navy', fontsize=11)
    plt.xlabel('cosine similarity')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, '+dataset+'+'+method+', Ngenes='+str(Ngenes))
    plt.savefig(fig_folder+'cos_sim/'+dataset_method+'_cos_sim_hist.png')
    plt.clf()
    # umapCompute
    adata_total_plot = adata_total.copy()
    adata_total_plot.obs['cos_sim'] = cos_sim
    scv.pl.velocity_embedding_stream(adata_total_plot, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                     save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_umapCompute.png")
    scv.pl.scatter(adata_total_plot, color='cos_sim', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_scatter_umapCompute.png")
    # umapOriginal
    adata_total_plot = adata_total.copy()
    adata_total_plot.obs['cos_sim'] = cos_sim
    adata_total_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    scv.pl.velocity_embedding_stream(adata_total_plot, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                     save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_umapOriginal.png")
    scv.pl.scatter(adata_total_plot, color='cos_sim', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_scatter_umapOriginal.png")

plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method)

######################################################
## plot velo_conf
scv.tl.velocity_confidence(total)
scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],
               save=fig_folder+"velo_conf/total_veloConf_umapCompute.png")

# plot velocity, velo_conf, cos_sim in one figure
## manually set vmin and vmax
total.obs['cos_sim'] = cos_sim
plt.clf()
fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))  # figsize=(horizontal, vertical)
scv.pl.velocity_embedding_stream(total, basis='umap',color="celltype",
                                 ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.2)
scv.pl.scatter(total,c='velocity_confidence',cmap='coolwarm',vmin=-.38,vmax=.97,ax=axs[1],legend_loc='none',frameon=False,size=100,alpha=0.15)
scv.pl.scatter(total,color='cos_sim',cmap='coolwarm',vmin=-.38,vmax=.97,ax=axs[2],legend_loc='none',frameon=False,size=100,alpha=0.15)
plt.savefig(fig_folder+"cos_sim/seed317_cosSim_veloConf.png")
plt.clf()

def plot_veloConf_and_cosSim(adata_total,adata_split1,adata_split2,adata_raw,dataset=dataset_short,method=method,fig_folder=fig_folder):
    if not 'velocity_confidence' in adata_total.obs.columns:
        scv.tl.velocity_confidence(adata_total)
    cos_sim = compute_cosine_similarity(adata_split1,adata_split2,method=method)[0]
    adata_total.obs['cos_sim'] = cos_sim
    vmin = np.min([0, np.min(cos_sim)-1e-5, np.min(adata_total.obs['velocity_confidence'])-1e-5])
    vmax = np.max([np.max(cos_sim)+1e-5, np.max(adata_total.obs['velocity_confidence'])+1e-5, 1])
    data_method = dataset+'_'+method
    # umapCompute
    scv.pl.scatter(adata_total, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"velo_conf/"+data_method+"_veloConf_umapCompute.png")
    plt.clf()
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                     ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.2)
    scv.pl.scatter(adata_total,c='velocity_confidence',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[1],legend_loc='none',frameon=False,size=100,alpha=0.15)
    scv.pl.scatter(adata_total,color='cos_sim',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[2],legend_loc='none',frameon=False,size=100,alpha=0.15)
    plt.savefig(fig_folder+"cos_sim/"+data_method+"_veloConf_and_cosSim_umapCompute.png")
    plt.clf()
    # umapOriginal
    adata_plot = adata_total.copy()
    adata_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    scv.pl.scatter(adata_plot, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"velo_conf/"+data_method+"_veloConf_umapOriginal.png")
    plt.clf()
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap',color="celltype",
                                     ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.2)
    scv.pl.scatter(adata_plot,c='velocity_confidence',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[1],legend_loc='none',frameon=False,size=100,alpha=0.15)
    scv.pl.scatter(adata_plot,color='cos_sim',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[2],legend_loc='none',frameon=False,size=100,alpha=0.15)
    plt.savefig(fig_folder+"cos_sim/"+data_method+"_veloConf_and_cosSim_umapOriginal.png")
    plt.clf()

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,method=method,fig_folder=fig_folder)

######################################################
## ptime
def plot_pseudotime(adata_in,data_raw,fig_name,dataset=dataset_short,method=method):
    data_method = dataset+'_'+method
    fig_name = data_method+'_'+fig_name
    if not 'velocity_confidence' in adata_in.obs.columns:
        scv.tl.velocity_confidence(adata_in)
    # umapCompute
    adata = adata_in.copy()
    scv.pl.scatter(adata, color="velocity_pseudotime", color_map="gnuplot", save=fig_folder+'ptime/'+fig_name+'_ptime_umapCompute.png')
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="celltype",ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.2)
    scv.pl.scatter(adata,ax=axs[1], color="velocity_pseudotime", color_map="gnuplot")
    plt.savefig(fig_folder+"ptime/"+fig_name+'_ptime_umapCompute_umap_and_ptime.png')
    plt.clf()
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = data_raw.obsm['X_umap'].copy()
    scv.pl.scatter(adata, color="velocity_pseudotime", color_map="gnuplot", save=fig_folder+'ptime/'+fig_name+'_ptime_umapOriginal.png')
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="celltype",ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.2)
    scv.pl.scatter(adata,ax=axs[1], color="velocity_pseudotime", color_map="gnuplot")
    plt.savefig(fig_folder+"ptime/"+fig_name+'_ptime_umapOriginal_umap_and_ptime.png')
    plt.clf()

plot_pseudotime(adata_in=split1,data_raw=raw,fig_name="split1",dataset=dataset_short,method=method)
plot_pseudotime(adata_in=split2,data_raw=raw,fig_name="split2",dataset=dataset_short,method=method)
plot_pseudotime(adata_in=total,data_raw=raw,fig_name="total",dataset=dataset_short,method=method)

def ptime_correlation_scatter_plot(s1,s2,seed,method,name,xlab,ylab,dataset="ery",cell_types=None,colors=None):
    if dataset=="ery":
        cell_types = total.obs['celltype']
        colors = dict(zip(['Blood progenitors 1', 'Blood progenitors 2', 'Erythroid1', 'Erythroid2', 'Erythroid3'],
                          ['#f9decf', '#c9a997', '#C72228', '#f79083', '#EF4E22']))
    data_method = dataset+'_'+method
    df = pd.DataFrame({'split1':s1.obs['velocity_pseudotime'], 'split2':s2.obs['velocity_pseudotime'],
                       'cell_types':cell_types})
    corr = np.round(np.corrcoef(s1.obs['velocity_pseudotime'],s2.obs['velocity_pseudotime']),3)[0,1]
    print(corr)
    plt.figure(figsize=(7, 5))
    for category, color in colors.items(): 
        print(category+' '+color)
        plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors))
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime correlation '+name+', corr='+str(corr)+')')
    plt.savefig(fig_folder+'ptime/'+data_method+"_pseudotimeCorr"+name+".png")
    plt.close()

ptime_correlation_scatter_plot(s1=split1,s2=split2,seed=317,method=method, name="split1vs2",xlab="split1",ylab="split2")
ptime_correlation_scatter_plot(s1=split1,s2=total,seed=317,method=method, name="split1vstotal",xlab="split1",ylab="total")
ptime_correlation_scatter_plot(s1=split2,s2=total,seed=317,method=method, name="split2vstotal",xlab="split2",ylab="total")


