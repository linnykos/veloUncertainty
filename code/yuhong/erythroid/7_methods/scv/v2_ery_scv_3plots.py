import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import datetime

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/scv/" 

total = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_total_v2.h5ad')
split1 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split1_v2.h5ad')
split2 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split2_v2.h5ad')

######################################################
## plot velocity - done in 1data.py, copied (and modified variable names here)
raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
scv.pl.velocity_embedding_stream(total, basis='umap',color="sequencing.batch",save=fig_folder+"velocity/ery_scv_total_umapCompute_byBatch.png")

def plot_velocities_scv(adata_in,adata_raw,fig_info,data,method,color_label=None):
    if data=="ery":
        color_label = 'celltype'
    elif data=="pan":
        color_label = 'clusters'
    data_method = data+"_"+method
    # umapCompute
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="celltype",save=fig_folder+"velocity/"+data_method+fig_info+"_umapCompute.png")
    # umapOriginal
    adata=adata_in.copy()
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=color_label,save=fig_folder+"velocity/"+data_method+fig_info+"_umapOriginal.png")    

plot_velocities_scv(total,raw,"total")
plot_velocities_scv(split1,raw,"seed317_split1")
plot_velocities_scv(split2,raw,"seed317_split2")

######################################################
## plot cosine similarity - done in 1data.py
velo_genes_split1 = split1.var.index[~np.isnan(split1.layers['velocity'][0])]
velo_genes_split2 = split2.var.index[~np.isnan(split2.layers['velocity'][0])]
common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0])) 

velo_df1 = pd.DataFrame(split1.layers['velocity'], columns=split1.var.index.tolist())
velo_df2 = pd.DataFrame(split2.layers['velocity'], columns=split2.var.index.tolist())
cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))
np.quantile(cos_sim,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

def compute_cosine_similarity_scv(adata_split1,adata_split2):
    velo_genes_split1 = adata_split1.var.index[~np.isnan(split1.layers['velocity'][0])]
    velo_genes_split2 = adata_split2.var.index[~np.isnan(split2.layers['velocity'][0])]
    common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0])) 
    velo_df1 = pd.DataFrame(adata_split1.layers['velocity'], columns=adata_split1.var.index.tolist())
    velo_df2 = pd.DataFrame(adata_split2.layers['velocity'], columns=adata_split2.var.index.tolist())
    cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))
    return cos_sim, common_genes_velocity.shape[0] # return cosine similarity and number of common genes in velocity computation

compute_cosine_similarity_scv(split1,split2)

# will be different
#cs1=np.diag(cosine_similarity(np.nan_to_num(adata_split1.layers['velocity'], nan=0),np.nan_to_num(adata_split2.layers['velocity'], nan=0)))
#np.quantile(cs1,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

import matplotlib.pyplot as plt
# cosine similarity
def plot_cosine_similarity_scv(adata_split1,adata_split2,adata_total, seed_split=317,text_x=None,text_y=None):
    cos_sim, Ngenes = compute_cosine_similarity_scv(adata_split1,adata_split2)
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
    plt.text(text_x, text_y, 'mean = '+str(np.round(np.mean(cos_sim),5)), color='navy', fontsize=11)
    plt.xlabel('cosine similarity (seed'+str(seed_split)+')')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, ery+scv, Ngenes='+str(Ngenes))
    plt.savefig(fig_folder+'cos_sim/seed'+str(seed_split)+'_cos_sim_hist.png')
    plt.clf()
    # umap
    adata_total_plot = adata_total.copy()
    adata_total_plot.obs['cos_sim'] = cos_sim
    scv.pl.velocity_embedding_stream(adata_total_plot, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                 save=fig_folder+"cos_sim/seed"+str(seed_split)+"_cos_sim_umapCompute.png")
    plt.clf()

plot_cosine_similarity_scv(split1,split2,total)

scv.tl.velocity_confidence(total)
scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],
               save=fig_folder+"velo_conf/total_veloConf_umapCompute.png")
scv.pl.scatter(total, color='cos_sim', cmap='coolwarm', perc=[1, 100],
               save=fig_folder+"cos_sim/seed317_cos_sim_umapCompute_scatter.png")

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


# ptime
scv.tl.velocity_pseudotime(split1)
scv.tl.velocity_pseudotime(split2)
scv.tl.velocity_pseudotime(total)
cell_types = total.obs['celltype']
colors = dict(zip(['Blood progenitors 1', 'Blood progenitors 2', 'Erythroid1', 'Erythroid2', 'Erythroid3'],
                  ['#f9decf', '#c9a997', '#C72228', '#f79083', '#EF4E22']))

def plot_pseudotime(adata_in,fig_name):
    adata = adata_in.copy()
    scv.tl.velocity_pseudotime(adata)
    scv.pl.scatter(adata, color="velocity_pseudotime", color_map="gnuplot", save=fig_folder+'ptime/'+fig_name+'_ptime_umapCompute.png')
    plt.clf()
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="celltype",ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.2)
    scv.pl.scatter(adata,ax=axs[1], color="velocity_pseudotime", color_map="gnuplot")
    plt.savefig(fig_folder+"ptime/"+fig_name+'_ptime_umapCompute_umap_and_ptime.png')
    plt.clf()

plot_pseudotime(split1,"seed317_split1")
plot_pseudotime(split2,"seed317_split2")
plot_pseudotime(total,"seed317_total")

def ptime_correlation_scatter_plot(s1,s2,seed,method,name,xlab,ylab,data="ery"):
    df = pd.DataFrame({'split1':s1.obs['velocity_pseudotime'], 'split2':s2.obs['velocity_pseudotime'],
                       'cell_types':cell_types})
    corr = np.round(np.corrcoef(s1.obs['velocity_pseudotime'],s2.obs['velocity_pseudotime']),3)[0,1]
    print(corr)
    data_method = data+"_"+method
    plt.figure(figsize=(7, 5))
    for category, color in colors.items(): 
        print(category+' '+color)
        plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors))
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime '+name+' (seed='+str(seed)+', corr='+str(corr)+')')
    plt.savefig(fig_folder+'ptime/'+data_method+"_seed"+str(seed)+"_pseudotimeCorr"+name+".png")
    plt.close()

ptime_correlation_scatter_plot(s1=split1,s2=split2,seed=317,method="scv", name="split1vs2",xlab="split1",ylab="split2")
ptime_correlation_scatter_plot(s1=split1,s2=total,seed=317,method="scv", name="split1vstotal",xlab="split1",ylab="total")
ptime_correlation_scatter_plot(s1=split2,s2=total,seed=317,method="scv", name="split2vstotal",xlab="split2",ylab="total")

