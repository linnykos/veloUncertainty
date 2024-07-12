import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'utv'
dataset_long = 'pancreas'
dataset_short = 'pan'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
raw = scv.read(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")

######################################################
## compute umap
#sc.tl.umap(total)
sc.tl.umap(split1)
sc.tl.umap(split2)

plot_velocities_scv_utv(adata_in=total,adata_raw=raw,fig_folder=fig_folder,fig_info="total",dataset=dataset_short,method=method)
plot_velocities_scv_utv(adata_in=split1,adata_raw=raw,fig_folder=fig_folder,fig_info="split1",dataset=dataset_short,method=method)
plot_velocities_scv_utv(adata_in=split2,adata_raw=raw,fig_folder=fig_folder,fig_info="split2",dataset=dataset_short,method=method)


######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)

plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## plot velo_conf
def plot_veloConf_and_cosSim_helper(adata_total,adata_raw,dataset,method,fig_folder,umapOriginal,vmin,vmax,Ngenes):
    adata_plot = adata_total.copy()
    celltype_label = None
    if dataset=="ery": celltype_label="celltype"
    elif dataset=="pan": celltype_label="clusters"
    data_method = dataset+'_'+method
    fig_umap = "umapCompute"
    if umapOriginal==True:
        adata_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
        fig_umap = "umapOriginal"
    scv.pl.scatter(adata_plot, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"velo_conf/"+data_method+"_veloConf_"+fig_umap+".png")
    plt.clf()
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata_plot,c='velocity_confidence',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[1],legend_loc='none',
                   title='velocity confidence, '+dataset+'+'+method,frameon=False,size=100,alpha=0.3)
    scv.pl.scatter(adata_plot,color='cos_sim',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[2],legend_loc='none',
                   title='velocity cosine similarity, '+dataset+'+'+method+', Ngenes='+str(Ngenes),frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+"cos_sim/"+data_method+"_veloConf_and_cosSim_"+fig_umap+".png")
    plt.clf()

def plot_veloConf_and_cosSim(adata_total,adata_split1,adata_split2,adata_raw,dataset=dataset_short,method=method,fig_folder=fig_folder):
    if not 'velocity_confidence' in adata_total.obs.columns:
        scv.tl.velocity_confidence(adata_total)
    celltype_label = None
    if dataset=="ery":
        celltype_label="celltype"
    elif dataset=="pan":
        celltype_label="clusters"
    cos_sim,Ngenes = compute_cosine_similarity(adata_split1,adata_split2,method=method)
    adata_total.obs['cos_sim'] = cos_sim
    vmin = np.min([0, np.min(cos_sim)-1e-5, np.min(adata_total.obs['velocity_confidence'])-1e-5])
    vmax = np.max([np.max(cos_sim)+1e-5, np.max(adata_total.obs['velocity_confidence'])+1e-5, 1])
    data_method = dataset+'_'+method
    # umapCompute
    plot_veloConf_and_cosSim_helper(adata_total,adata_raw,dataset,method,fig_folder,umapOriginal=False,vmin=vmin,vmax=vmax,Ngenes=Ngenes)
    # umapOriginal
    plot_veloConf_and_cosSim_helper(adata_total,adata_raw,dataset,method,fig_folder,umapOriginal=True,vmin=vmin,vmax=vmax,Ngenes=Ngenes)

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## ptime
'velocity_pseudotime' in split1.obs.columns
#scv.tl.velocity_pseudotime(total)
#scv.tl.velocity_pseudotime(split1)
#scv.tl.velocity_pseudotime(split2)

plot_pseudotime(adata_in=split1,data_raw=raw,fig_name="split1",dataset=dataset_short,method=method)
plot_pseudotime(adata_in=split2,data_raw=raw,fig_name="split2",dataset=dataset_short,method=method)
plot_pseudotime(adata_in=total,data_raw=raw,fig_name="total",dataset=dataset_short,method=method)

ptime_correlation_scatter_plot(s1=split1,s2=split2,method=method,dataset='pan',name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split1,s2=total,method=method,dataset='pan',name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split2,s2=total,method=method,dataset='pan',name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder)


