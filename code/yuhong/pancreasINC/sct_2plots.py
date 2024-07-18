import sctour as sct
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from sctour_misc import *
from v2_functions import *

method = 'sct'
dataset_long = 'pancreasINC'
dataset_short = 'panINC'

sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_pancreasINC/sct/adata_panINC_sct_total_v2.h5ad')
split1 = sc.read_h5ad(data_folder+'v2_pancreasINC/sct/adata_panINC_sct_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_pancreasINC/sct/adata_panINC_sct_split2_v2.h5ad')

raw = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
cell_index = np.array(np.where(raw.obs['clusters']!='Pre-endocrine')[0])
def create_adata_INC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}#adata_old.uns['clusters_colors'].copy()
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

S_raw = raw.layers['spliced'][cell_index,:]
U_raw = raw.layers['unspliced'][cell_index,:]
raw = create_adata_INC(S=S_raw,U=U_raw,adata_old=raw)

total.uns['clusters_colors'] = split1.uns['clusters_colors'].copy()
######################################################
#### plot scTour vector field
def plot_vf_umap(adata_in,adata_raw,fig_name,data,method='sct'):
    celltype_label = "celltype"
    if "pan" in data: celltype_label="clusters"
    data_method = data+'_'+method
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color=celltype_label, ax=axs[0], legend_loc='on data', show=False, frameon=False)
    print("umapOriginal[0] done")
    sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
    print("umapOriginal[1] done")
    sct.vf.plot_vector_field(adata,zs_key='X_TNODE',vf_key='X_VF',use_rep_neigh='X_TNODE',color=celltype_label, 
                             show=False,ax=axs[2],legend_loc='none',frameon=False,size=100,alpha=0.2,title=data+' '+fig_name,
                             save=fig_folder+'vf/'+data_method+'_vf_'+fig_name+'_umapOriginal.png')
    print("umapOriginal[2] done")    
    # umapCompute
    adata = adata_in.copy()
    adata = adata[np.argsort(adata.obs['ptime'].values), :]
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=15) # used 30 in first ery version
    sc.tl.umap(adata)
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color=celltype_label, ax=axs[0], legend_loc='on data', show=False, frameon=False)
    print("umapCompute[0] done")
    sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
    print("umapCompute[1] done")
    sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color=celltype_label, 
                            show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, title=data+' '+fig_name,
                            save=fig_folder+'vf/'+data_method+'_vf_'+fig_name+'_umapCompute.png')
    print("umapCompute[2] done")
 
plot_vf_umap(adata_in=split1, adata_raw=raw, fig_name="split1",data=dataset_short,method=method)
plot_vf_umap(adata_in=split2, adata_raw=raw, fig_name="split2",data=dataset_short,method=method)
plot_vf_umap(adata_in=total, adata_raw=raw, fig_name="total",data=dataset_short,method=method)

######################################################
#### plot velocity
def plot_velocity_sct(adata_in, adata_raw, fig_name, dataset, fig_folder, method='sct'):
    data_method = dataset+'_'+method
    print(data_method)
    celltype_label = "celltype"
    if "pan" in dataset: celltype_label = 'clusters'
    if not fig_name == "total":
        # umapOriginal: for total: ValueError: Your neighbor graph seems to be corrupted. Consider recomputing via pp.neighbors.
        adata = adata_in.copy()
        adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
        scv.tl.velocity_graph(adata)
        scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label, title='Velocity '+dataset+'+'+method+' '+fig_name,
                                        save=fig_folder+"velocity/"+data_method+"_"+fig_name+"_umapOriginal.png")
    # umapOriginal_recomputeNbr
    adata = adata_in.copy()
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata)    
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=10) # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label, title='Velocity '+dataset+'+'+method+' '+fig_name,
                                     save=fig_folder+"velocity/"+data_method+"_"+fig_name+"_umapOriginal_recomputeNbr.png")
    # umapCompute
    adata = adata_in.copy()
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata)    
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=10) # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label, title='Velocity '+dataset+'+'+method+' '+fig_name,
                                     save=fig_folder+"velocity/"+data_method+"_"+fig_name+"_umapCompute.png")


plot_velocity_sct(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,fig_folder=fig_folder)
plot_velocity_sct(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,fig_folder=fig_folder)
plot_velocity_sct(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,fig_folder=fig_folder)


######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,text_x=0)

plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)


######################################################
## plot velo_conf
plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_veloConf_hist(adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## ptime
plot_pseudotime(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,method=method,fig_folder=fig_folder)

ptime_correlation_scatter_plot(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder)


