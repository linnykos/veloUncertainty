import sctour as sct
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime
import os
from sklearn.metrics.pairwise import cosine_similarity

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/" 

### read in data
split1 = sc.read_h5ad(data_folder+'v2_erythroid/sct/adata_ery_sct_seed317_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_erythroid/sct/adata_ery_sct_seed317_split2_v2.h5ad')
total = sc.read_h5ad(data_folder+'v2_erythroid/sct/adata_ery_sct_total_v2.h5ad')
raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")


#### plot vf
save_path = os.path.join(save_dir, "Writeup6_eprs_mg_scvi_ePRS.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    adata2,
    color=["ePRS"],
    frameon=False,
    title="By ePRS status",
    size=5,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

def plot_vf_umap(adata_in,adata_raw,fig_name):
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    #save_path_ptime = os.path.join(fig_folder, "sct/ptime/ery_sct_ptime_"+fig_name+"_umapOriginal.png")
    sc.pl.umap(adata, color='celltype', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
    sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', 
                            show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                            save=fig_folder+'sct/vf/ery_sct_vf_'+fig_name+'_umapOriginal.png')
    adata = adata_in.copy()
    adata = adata[np.argsort(adata.obs['ptime'].values), :]
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(adata)
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color='celltype', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
    sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', 
                            show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                            save=fig_folder+'sct/vf/ery_sct_vf_'+fig_name+'_umapCompute.png')
 
plot_vf_umap(adata_in=split1, adata_raw=raw, fig_name="seed317_split1")
plot_vf_umap(adata_in=split2, adata_raw=raw, fig_name="seed317_split2")
plot_vf_umap(adata_in=total, adata_raw=raw, fig_name="total")

#### plot velocity
def plot_velocity(adata_in, adata_raw, fig_name):
    # umapOriginal
    data = adata_in.copy()
    data.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis='umap',color="celltype",save=fig_folder+"sct/velocity/ery_sct_"+fig_name+"_umapOriginal.png")
    # umapOriginal_recomputeNbr
    data = adata_in.copy()
    data.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    sc.pp.neighbors(data, use_rep='X_TNODE', n_neighbors=30)
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis='umap',color="celltype",save=fig_folder+"sct/velocity/ery_sct_"+fig_name+"_umapOriginal_recomputeNbr.png")
    # umapCompute
    data = adata_in.copy()
    sc.pp.neighbors(data, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(data)
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis='umap',color="celltype",save=fig_folder+"sct/velocity/ery_sct_"+fig_name+"_umapCompute.png")

plot_velocity(adata_in=split1,adata_raw=raw, fig_name="seed317_split1")
plot_velocity(adata_in=split2, adata_raw=raw, fig_name="seed317_split2")
plot_velocity(adata_in=total, adata_raw=raw, fig_name="total")


#### plot cosine similarity
### common genes being selected and used for velocity computation
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index))
common_genes.shape # 1304->1237 common genes
np.intersect1d(np.array(split1.var.index), np.array(total.var.index)).shape # 229 -> 1531
np.intersect1d(np.array(split2.var.index), np.array(total.var.index)).shape # 244 -> 1519

velo_df1 = pd.DataFrame(split1.layers['velocity'], columns=split1.var.index.tolist())
velo_df2 = pd.DataFrame(split2.layers['velocity'], columns=split2.var.index.tolist())
cos_sim = np.diag(cosine_similarity(velo_df1[common_genes],velo_df2[common_genes]))
np.quantile(cos_sim,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

def plot_cosine_similarity_sct(cos_sim, total,adata_raw, seed_split=317,text_x=None,text_y=None):
    # histogram
    plt.clf()
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
    plt.title('Histogram of cosine similarity, ery+sct')
    plt.savefig(fig_folder+'sct/cos_sim/seed'+str(seed_split)+'_cos_sim_hist.png')
    plt.clf()
    # umap
    adata = total.copy()
    adata.obs['cos_sim'] = cos_sim
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                 save=fig_folder+"sct/cos_sim/seed"+str(seed_split)+"_cos_sim_umapCompute.png")
    scv.pl.scatter(adata, color='cos_sim', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"sct/cos_sim/seed"+str(seed_split)+"_cos_sim_umapCompute_scatter.png")
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color='celltype', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                     ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2,
                                     save=fig_folder+"sct/cos_sim/seed"+str(seed_split)+"_cos_sim_umapCompute_wCelltype.png")
    
    adata = total.copy()
    adata.obs['cos_sim'] = cos_sim
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                 save=fig_folder+"sct/cos_sim/seed"+str(seed_split)+"_cos_sim_umapOriginal.png")
    scv.pl.scatter(adata, color='cos_sim', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"sct/cos_sim/seed"+str(seed_split)+"_cos_sim_umapOriginal_scatter.png")


plot_cosine_similarity_sct(cos_sim,total,adata_raw=raw)

def plot_cosine_similarity_sct_wCelltype(cos_sim, total,adata_raw, seed_split=317,text_x=None,text_y=None):
    # umap
    adata = total.copy()
    adata.obs['cos_sim'] = cos_sim
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color='celltype', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                     ax=axs[1], legend_loc='none', frameon=False, size=100, alpha=0.2,
                                     save=fig_folder+"sct/cos_sim/seed"+str(seed_split)+"_cos_sim_umapCompute_wCelltype.png")
    adata = total.copy()
    adata.obs['cos_sim'] = cos_sim
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    scv.tl.velocity_graph(adata)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color='celltype', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                     ax=axs[1], legend_loc='none', frameon=False, size=100, alpha=0.2,
                                     save=fig_folder+"sct/cos_sim/seed"+str(seed_split)+"_cos_sim_umapOriginal_wCelltype.png")

plot_cosine_similarity_sct_wCelltype(cos_sim,total,adata_raw=raw)



#### plot veloConf
scv.tl.velocity_graph(total)
scv.tl.velocity_confidence(total)
# KeyError: 'Ms'
scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],
               save=fig_folder+"sct/velo_conf/total_veloConf_umapOriginal.png")

#### plot ptime corr
cell_to_color = dict(zip(np.array(raw.obs['celltype'].cat.categories), raw.uns['celltype_colors']))
cell_types = raw.obs['celltype']

def plot_ptime_corr_scatter(adata_split1,adata_split2,seed_split,name,xlab,ylab):
    s1 = adata_split1.copy()
    s2 = adata_split2.copy()
    df = pd.DataFrame({'split1':np.array(s1.obs['ptime']), 'split2':np.array(s2.obs['ptime']),
                       'cell_types':np.array(cell_types)})
    corr = np.round(np.corrcoef(s1.obs['ptime'],s2.obs['ptime']),3)[0,1]
    for category, color in cell_to_color.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(cell_to_color))
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime '+name+' ery+sct (seed='+str(seed_split)+', corr='+str(corr)+')')
    plt.savefig(fig_folder+"sct/ptime/seed"+str(seed_split)+"_pseudotime"+name+".png")
    plt.close()

plot_ptime_corr_scatter(adata_split1=split1,adata_split2=split2,seed_split=317,name="split1vs2",xlab="split1",ylab="split2")
plot_ptime_corr_scatter(adata_split1=split1,adata_split2=total,seed_split=317,name="split1vstotal",xlab="split1",ylab="total")
plot_ptime_corr_scatter(adata_split1=split2,adata_split2=total,seed_split=317,name="split2vstotal",xlab="split2",ylab="total")




