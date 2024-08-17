import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import scvelo as scv
import datetime

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v1_erythroid/sct/" 

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from sctour_misc import *

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

#########################################################
# read in data
adata = sc.read_h5ad(data_folder+"v1_erythroid/sct/adata_ery_sct_preprocess.h5ad")
tnode = sct.predict.load_model(data_folder+'v1_erythroid/sct/tnode_ery_sct_preprocess.pth')

adata_split1 = sc.read_h5ad(data_folder+'v1_erythroid/sct/adata_ery_sct_seed317_split1.h5ad')
tnode_split1 = sct.predict.load_model(data_folder+'v1_erythroid/sct/tnode_ery_sct_seed317_split1.pth')

adata_split2 = sc.read_h5ad(data_folder+'v1_erythroid/sct/adata_ery_sct_seed317_split2.h5ad')
tnode_split2 = sct.predict.load_model(data_folder+'v1_erythroid/sct/tnode_ery_sct_seed317_split2.pth')



#########################################################
# plot velocities
def plot_velocity(split, total, fig_info):
    print_message_with_time("##################### Calling plot function")
    data = split.copy()
    #data = sc.read_h5ad(data_folder+'v1_erythroid/sct/adata_seed'+str(seed_split)+'_split1.h5ad')
    data.obsm['X_umap'] = total.obsm['X_umap'].copy()
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis='umap',color="celltype",
                                     save=fig_folder+"velocity/ery_sct_"+fig_info+"_umapOriginal.png")
    print_message_with_time("##################### umapOriginal plotted")
    data = split.copy()
    #data = sc.read_h5ad(data_folder+'v1_erythroid/sct/adata_seed'+str(seed_split)+'_split1.h5ad')
    sc.pp.neighbors(data, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(data)
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis='umap',color="celltype",
                                     save=fig_folder+"velocity/ery_sct_"+fig_info+"_umapCompute.png")
    print_message_with_time("##################### umapCompute plotted")
    

plot_velocity(split=adata_split1,total=adata,fig_info="seed317_split1")
plot_velocity(split=adata_split2,total=adata,fig_info="seed317_split2")

sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
# do this because of the error:
# ValueError: Your neighbor graph seems to be corrupted. Consider recomputing via pp.neighbors.
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap',color="clusters", save=fig_folder+"velocity/ery_sct_preprocess.png")

#########################################################
# plot scTour vector field
def plot_vf_umap(res,name):
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(16, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(res, color='celltype', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    sc.pl.umap(res, color='ptime', ax=axs[1], show=False, frameon=False)
    sct.vf.plot_vector_field(res, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', 
                            show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                            save=fig_folder+"vf/"+name)

adata_split1.obsm['X_umap'] = adata.obsm['X_umap'].copy()
adata_split2.obsm['X_umap'] = adata.obsm['X_umap'].copy()

plot_vf_umap(res=adata_split1, name="seed317_split1_umapOriginal.png")
plot_vf_umap(res=adata_split2, name="seed317_split2_umapOriginal.png")

def compute_umap_and_plot_vf(res,name):
    del res.obsm['X_umap']
    res = res[np.argsort(res.obs['ptime'].values), :]
    sc.pp.neighbors(res, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(res)
    plot_vf_umap(res=res, name=name)

compute_umap_and_plot_vf(res=adata_split1, name="seed317_split1_umapCompute.png")
compute_umap_and_plot_vf(res=adata_split2, name="seed317_split2_umapCompute.png")


#########################################################
# cosine similarity
def plot_cosine_similarity(split1, split2, total, seed_split,text_x=None,text_y=None):
    cos_sim = np.diag(cosine_similarity(split1.layers['velocity'],split2.layers['velocity']))
    # histogram
    plt.clf()
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='black') 
    max_frequency = np.max(counts)
    if text_x is None:
        text_x = np.quantile(cos_sim,[.05])
    if text_y is None:
        text_y = max_frequency/2
    plt.axvline(np.mean(cos_sim), color='red', linestyle='dashed', linewidth=1) ## add mean
    plt.text(text_x, text_y, 'mean = '+str(np.round(np.mean(cos_sim),4)), color='blue', fontsize=10)
    plt.xlabel('cosine similarity (seed'+str(seed_split)+')')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, ery+sct')
    plt.savefig(fig_folder+'cos_sim/seed'+str(seed_split)+'_cos_sim_hist.png')
    plt.clf()
    # umap
    total.obs['cos_sim_seed'] = pd.DataFrame(cos_sim)
    sct.vf.plot_vector_field(total, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='cos_sim_seed'+str(seed_split), 
                             cmap='coolwarm',legend_loc='none', frameon=False, size=100, alpha=0.2, 
                             save=fig_folder+"cos_sim/seed"+str(seed_split)+"_cos_sim_vf.png")
    scv.pl.velocity_embedding_stream(total, basis='umap',color='cos_sim_seed'+str(seed_split),
                                     save=fig_folder+"velocity/seed"+str(seed_split)+"_umap.png")
    scv.pl.velocity_embedding_stream(total, basis='pca',color='cos_sim_seed'+str(seed_split),
                                     save=fig_folder+"velocity/seed"+str(seed_split)+"_pca.png")
    del total.obs['cos_sim_seed']


#########################################################
### velo confidence

#########################################################
### pseudotime
## pseudotime in these objects have been reordered for the purpose of plotting

cell_to_color = dict(zip(np.array(adata.obs['celltype'].cat.categories), adata.uns['celltype_colors']))
cell_types = adata.obs['celltype']

def ptime_scatter_plot(s1,s2,seed_split,name,xlab,ylab):
    df = pd.DataFrame({'split1':np.array(s1.obs['ptime']), 'split2':np.array(s2.obs['ptime']),
                       'cell_types':np.array(cell_types)})
    corr = np.round(np.corrcoef(s1.obs['ptime'],s2.obs['ptime']),3)[0,1]
    for category, color in cell_to_color.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(cell_to_color))
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime '+name+' (seed='+str(seed_split)+', corr='+str(corr)+')')
    plt.savefig(fig_folder+"ptime/"+"seed"+str(seed_split)+"_pseudotime"+name+".png")
    plt.close()

ptime_scatter_plot(s1=adata_split1,s2=adata_split2,seed_split=317,name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=adata_split1,s2=adata,seed_split=317,name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=adata_split2,s2=adata,seed_split=317,name="split2vstotal",xlab="split2",ylab="total")


