import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import random
import matplotlib.pyplot as plt
import sctour as sct

sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from sctour_misc import *
from v4_functions import print_message_with_time,read_raw_adata,read_data_v4,compute_cosine_similarity_union,get_celltype_label,get_basis_type,get_metric_color_and_title,read_raw_adata

import datetime
def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

############################################
### train sct model
def sct_train_and_return_tnode(adata, sct_seed=615):
    import sctour as sct
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    adata.X = adata.X.astype(np.float32)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)
    tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
    return tnode

from scipy.stats import spearmanr
def test_timestep(adata_split1,adata_split2,adata_total,tnode1,tnode2,tnode,time,sct_seed=615):
    #split1 = adata_split1.copy()
    #split2 = adata_split2.copy()
    total = adata_total.copy()
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    #split1.layers['velocity'] = compute_sctour_velocity(tnode1, timestep=time)
    #split2.layers['velocity'] = compute_sctour_velocity(tnode2, timestep=time)
    total.layers['velocity'] = compute_sctour_velocity(tnode, timestep=time)
    #cos_sim,Ngenes = compute_cosine_similarity_union(split1,split2,method='sct')
    scv.tl.velocity_graph(total,n_jobs=8)
    scv.tl.velocity_pseudotime(total)
    #scv.tl.recover_dynamics(total,var_names=total.var.index,n_jobs=8)
    #scv.tl.latent_time(total)
    ptime_cor = spearmanr(total.obs['ptime'], total.obs['velocity_pseudotime']).correlation
    #latent_cor = spearmanr(total.obs['ptime'], total.obs['latent_time']).correlation
    #print([np.mean(cos_sim), np.median(cos_sim)])
    print(ptime_cor)
    #print([ptime_cor, latent_cor])
    #return np.round(np.mean(cos_sim),6),np.round(np.median(cos_sim),6),ptime_cor,latent_cor
    #return np.round(np.mean(cos_sim),6),np.round(np.median(cos_sim),6),ptime_cor
    return np.round(np.mean(ptime_cor),5)

def read_data_and_run_sct(dataset_short,dataset_long,method,data_folder,savedata_folder,split_version,split_seed,sct_seed):
    print_message_with_time("################## Read data")
    raw = read_raw_adata(dataset_short)
    adata = read_data_v4(dataset_long,dataset_short,method,split_seed,split_version,allgenes=True,outputAdded=False)
    gene_names = adata.var.index.copy()
    positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
    #S_mat = adata.layers['spliced'].copy()
    #U_mat = adata.layers['unspliced'].copy()
    if (not 'larry' in dataset_long) and (not dataset_long=='greenleaf'):
        adata.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
        adata.obsm['X_pcaOriginal'] = raw.obsm['X_pca'].copy()
    print_message_with_time("########### start to train model for "+split_version+' ')
    tnode = sct_train_and_return_tnode(adata, sct_seed)
    print_message_with_time("########### start to compute velocity for "+split_version+' ')
    print_message_with_time("########### "+split_version+' velocity computed, start to write data')
    adata.layers['spliced_original'] = adata.layers['spliced'].copy()
    adata.layers['unspliced_original'] = adata.layers['unspliced'].copy()
    adata.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+split_version+'_v4.h5ad')
    tnode.save_model(save_dir=savedata_folder, save_prefix='tnode_'+dataset_short+'_'+method+'_'+split_version+'_v4')
    print_message_with_time("########### "+split_version+' data wrote')

def sct_add_state_info_colors(adata,colors=["#6e8ea1","#ffab6e","#91b68a","#d18887","#a892c4","#b2967e","#dba8bc","#a0a0a0","#c4c88a","#87c3c9","#ffcccc"]):
    adata.uns['state_info_colors'] = colors

########################################
def get_umap_sct(adata,umapOriginal=False,moments=True,velocity_graph=True):
    if moments==True:
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
    if umapOriginal==True:
        adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    else:
        sc.tl.umap(adata)
    if velocity_graph==True: 
        scv.tl.velocity_graph(adata,n_jobs=8)

def plot_vf_umap(adata_in,data_version,data,fig_folder,method='sct',celltype_label=None):
    if celltype_label==None: celltype_label = get_celltype_label(data)
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    #scv.tl.velocity_graph(adata)
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color=celltype_label, ax=axs[0], legend_loc='on data', show=False, frameon=False)
    print("umapOriginal[0] done")
    sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
    print("umapOriginal[1] done")
    sct.vf.plot_vector_field(adata,zs_key='X_TNODE',vf_key='X_VF',use_rep_neigh='X_TNODE',color=celltype_label, 
                             show=False,ax=axs[2],legend_loc='none',frameon=False,size=100,alpha=0.2,title=data+' '+data_version,
                             save=fig_folder+'vf/'+data+'_'+method+'_vf_'+data_version+'_umapOriginal.png')
    print("umapOriginal[2] done")    
    # umapCompute
    adata = adata_in.copy()
    adata = adata[np.argsort(adata.obs['ptime'].values), :]
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30) # used 30 in first ery version
    sc.tl.umap(adata)
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color=celltype_label, ax=axs[0], legend_loc='on data', show=False, frameon=False)
    print("umapCompute[0] done")
    sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
    print("umapCompute[1] done")
    sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color=celltype_label, 
                            show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, title=data+' '+data_version,
                            save=fig_folder+'vf/'+data+'_'+method+'_vf_'+data_version+'_umapCompute.png')
    print("umapCompute[2] done")
 

def plot_sct_velocity(adata_in,data_version,dataset,fig_folder,recompute=True,method='sct',celltype_label=None):
    data_method = dataset+'_'+method
    print(data_method)
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    scv.tl.velocity_graph(adata,n_jobs=8)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,recompute=recompute, 
                                        title='Velocity '+dataset+'+'+method+' '+data_version,
                                        save=fig_folder+"velocity/"+data_method+"_"+data_version+"_umapOriginal.png")
    # umapCompute
    adata = adata_in.copy()
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,recompute=recompute,
                                     title='Velocity '+dataset+'+'+method+' '+data_version,
                                     save=fig_folder+"velocity/"+data_method+"_"+data_version+"_umapCompute.png")


