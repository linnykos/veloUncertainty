# run this in sctClone
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import sctour as sct

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'sct'
dataset_long = 'pancreas'
dataset_short = 'pan'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
raw = scv.read(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")

######################################################
#### plot scTour vector field
def plot_vf_umap(adata_in,adata_raw,fig_name,data,method='sct'):
    celltype_label = "celltype"
    if data=="pan": celltype_label="clusters"
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
plot_velocity_sct(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,fig_folder=fig_folder,recompute=True)
plot_velocity_sct(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,fig_folder=fig_folder,recompute=True)
plot_velocity_sct(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,fig_folder=fig_folder,recompute=True)

plot_velocity_sct(adata_in=split1,adata_raw=raw,fig_name="recompF_split1",dataset=dataset_short,fig_folder=fig_folder,recompute=False)
plot_velocity_sct(adata_in=split2,adata_raw=raw,fig_name="recompF_split2",dataset=dataset_short,fig_folder=fig_folder,recompute=False)
plot_velocity_sct(adata_in=total,adata_raw=raw,fig_name="recompF_total",dataset=dataset_short,fig_folder=fig_folder,recompute=False)


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

ptime_correlation_scatter_plot(s1=split1,s2=split2,method=method,dataset='pan',name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split1,s2=total,method=method,dataset='pan',name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split2,s2=total,method=method,dataset='pan',name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder)



