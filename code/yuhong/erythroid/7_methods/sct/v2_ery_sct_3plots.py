# run this in sctClone
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import sctour as sct

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'sct'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 9815 × 2000
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split2_v2.h5ad') # 
raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad") # 9815 × 53801

######################################################
#### plot vf
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
plot_velocity_sct(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,fig_folder=fig_folder)
plot_velocity_sct(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,fig_folder=fig_folder)
plot_velocity_sct(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,fig_folder=fig_folder)


######################################################
## plot cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,text_x=0)

plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)



######################################################
#### plot veloConf
#scv.tl.velocity_graph(total)
#scv.tl.velocity_confidence(total)
# KeyError: 'Ms'
#scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],save=fig_folder+"sct/velo_conf/total_veloConf_umapOriginal.png")

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)

######################################################
## ptime
plot_pseudotime(adata_in=split1,adata_raw=raw,fig_name="split1",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=split2,adata_raw=raw,fig_name="split2",dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_pseudotime(adata_in=total,adata_raw=raw,fig_name="total",dataset=dataset_short,method=method,fig_folder=fig_folder)

ptime_correlation_scatter_plot(s1=split1,s2=split2,method=method,dataset=dataset_short,name="split1vs2",xlab="split1",ylab="split2",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split1,s2=total,method=method,dataset=dataset_short,name="split1vstotal",xlab="split1",ylab="total",fig_folder=fig_folder)
ptime_correlation_scatter_plot(s1=split2,s2=total,method=method,dataset=dataset_short,name="split2vstotal",xlab="split2",ylab="total",fig_folder=fig_folder)


exit()

######################################################
######################################################
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

######################################################
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




