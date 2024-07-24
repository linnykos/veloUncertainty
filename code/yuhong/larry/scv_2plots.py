import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
import anndata as ad
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'scv'
dataset_long = 'larry'
dataset_short = 'larry'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
raw = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/larry.h5ad') # n_obs × n_vars = 49302 × 23420
raw.obsm['X_umap']=raw.obsm['X_emb']

import matplotlib as mpl
def rgb2hex(rgb):
    r = int(rgb[0]*255)
    g = int(rgb[1]*255)
    b = int(rgb[2]*255)
    return "#{:02x}{:02x}{:02x}".format(r,g,b)

split1.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]
split2.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]
total.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]

######################################################
## plot velocity
plot_velocity_scv_utv(adata_in=total,adata_raw=raw,fig_folder=fig_folder,fig_info="total",dataset=dataset_short,method=method)
plot_velocity_scv_utv(adata_in=split1,adata_raw=raw,fig_folder=fig_folder,fig_info="split1",dataset=dataset_short,method=method)
plot_velocity_scv_utv(adata_in=split2,adata_raw=raw,fig_folder=fig_folder,fig_info="split2",dataset=dataset_short,method=method)

plot_velocity_scv_utv(adata_in=total,adata_raw=raw,fig_folder=fig_folder,fig_info="total_recompF",dataset=dataset_short,method=method,recompute=False)

######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
""" Number of overlapped genes for velocity computation in splits = 624 """
plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_hist_by_celltype(split1,split2,total,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## plot velo_conf
plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_veloConf_hist(adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## ptime
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total)
    scv.tl.velocity_pseudotime(split1)
    scv.tl.velocity_pseudotime(split2)
"""
WARNING: Uncertain or fuzzy root cell identification. Please verify.
    identified 10 regions of root cells and 3 regions of end points .
"""
plot_pseudotime(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,method=method,fig_folder=fig_folder)

ptime_correlation_scatter_plot(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder)


