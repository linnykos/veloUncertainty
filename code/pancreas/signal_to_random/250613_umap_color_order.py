import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
import os
import re

split_seed = 317


import matplotlib.pyplot as plt
import numpy as np

def plot_velo_conf_cos_sim_matpltlib_sorted(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,celltype_label=None, rev=False):
    if (not 'velocity_confidence' in adata_total.obs.columns): scv.tl.velocity_confidence(adata_total)
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    vmin = np.min( [np.min( adata_total.obs['velocity_confidence'] ), np.min( adata_total.obs['cos_sim'] )] )
    vmax = np.max( [np.max( adata_total.obs['velocity_confidence'] ), np.max( adata_total.obs['cos_sim'] )] )
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    # velo_conf
    adata = adata_total.copy()
    X = adata.obsm["X_umapOriginal"][:, 0]
    Y = adata.obsm["X_umapOriginal"][:, 1]
    color_values = adata_total.obs['velocity_confidence']
    if rev: sort_idx = np.argsort(color_values)[::-1]
    else: sort_idx = np.argsort(color_values)  # Smallest at the end
    X_sorted = X[sort_idx]
    Y_sorted = Y[sort_idx]
    color_sorted = color_values[sort_idx]
    # Plot: velo_conf
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sc = ax.scatter( X_sorted, Y_sorted, c=color_sorted, cmap="coolwarm", vmin=vmin, vmax=vmax, s=6, alpha=0.6, edgecolors="none" )
    ax.set_title("Local coherence " + dataset + '+' + method + ' (Ngenes=' + str(adata.shape[1]) + ')')
    ax.axis("off")
    plt.colorbar(sc, ax=ax, label='local coherence')
    plt.tight_layout()
    if rev: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_velo_conf_umapOriginal_sorted_matpltlib-1.png')
    else: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_velo_conf_umapOriginal_sorted_matpltlib.png')
    plt.clf()
    # 
    color_values = cos_sim
    if rev: sort_idx = np.argsort(color_values)[::-1]
    else: sort_idx = np.argsort(color_values)  # Smallest at the end
    X_sorted = X[sort_idx]
    Y_sorted = Y[sort_idx]
    color_sorted = color_values[sort_idx]
    # Plot: cos_sim
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sc = ax.scatter( X_sorted, Y_sorted, c=color_sorted, cmap="coolwarm", vmin=vmin, vmax=vmax, s=6, alpha=0.6, edgecolors="none" )
    ax.set_title("Replicate coherence " + dataset + '+' + method + ' (Ngenes=' + str(Ngenes) + ')')
    ax.axis("off")
    plt.colorbar(sc, ax=ax, label='cos_sim')
    plt.tight_layout()
    if rev: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_cos_sim_umapOriginal_sorted_matpltlib-1.png')
    else: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_cos_sim_umapOriginal_sorted_matpltlib.png')
    plt.clf()


def plot_velo_conf_matpltlib_sorted(adata_total,dataset,method,fig_folder,split_seed,celltype_label=None,rev=False):
    if (not 'velocity_confidence' in adata_total.obs.columns): scv.tl.velocity_confidence(adata_total)
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    adata = adata_total.copy()
    X = adata.obsm["X_umapOriginal"][:, 0]
    Y = adata.obsm["X_umapOriginal"][:, 1]
    color_values = adata_total.obs['velocity_confidence']
    if rev: sort_idx = np.argsort(color_values)[::-1]
    else: sort_idx = np.argsort(color_values)  # Smallest at the end
    X_sorted = X[sort_idx]
    Y_sorted = Y[sort_idx]
    color_sorted = color_values[sort_idx]
    # Plot
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sc = ax.scatter( X_sorted, Y_sorted, c=color_sorted, cmap="coolwarm", vmin=-1, vmax=1, s=5, alpha=0.5, edgecolors="none" )
    ax.set_title("Local coherence " + dataset + '+' + method + ' (Ngenes=' + str(adata.shape[1]) + ')')
    ax.axis("off")
    plt.colorbar(sc, ax=ax, label='local coherence')
    plt.tight_layout()
    if rev: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_velo_conf_umapOriginal_sorted_matpltlib-1.png')
    else: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_velo_conf_umapOriginal_sorted_matpltlib.png')
    plt.clf()

def plot_cosine_similarity_withRef_matpltlib_sorted(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,celltype_label=None, rev=False):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    adata = adata_total.copy()
    X = adata.obsm["X_umapOriginal"][:, 0]
    Y = adata.obsm["X_umapOriginal"][:, 1]
    color_values = cos_sim
    if rev: sort_idx = np.argsort(color_values)[::-1]
    else: sort_idx = np.argsort(color_values)  # Smallest at the end
    X_sorted = X[sort_idx]
    Y_sorted = Y[sort_idx]
    color_sorted = color_values[sort_idx]
    # Plot
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sc = ax.scatter( X_sorted, Y_sorted, c=color_sorted, cmap="coolwarm", vmin=-1, vmax=1, s=5, alpha=0.5, edgecolors="none" )
    ax.set_title("Replicate coherence " + dataset + '+' + method + ' (Ngenes=' + str(Ngenes) + ')')
    ax.axis("off")
    plt.colorbar(sc, ax=ax, label='cos_sim')
    plt.tight_layout()
    if rev: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_cos_sim_umapOriginal_sorted_matpltlib-1.png')
    else: plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_cos_sim_umapOriginal_sorted_matpltlib.png')
    plt.clf()


def plot_umap_color_sorted_pancreas(method, split_seed=317, dataset_short='pan', dataset_long = 'pancreas'):
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method+'/')
    if 'scv' in method or 'utv' in method:
        path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1', s)][0]
        path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2', s)][0]
        path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total', s)][0]
    else: 
        path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1.*outputAdded', s)][0]
        path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2.*outputAdded', s)][0]
        path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total.*outputAdded', s)][0]
    split1 = sc.read_h5ad( path_split1 )
    split2 = sc.read_h5ad( path_split2 )
    total = sc.read_h5ad( path_total )
    plot_velo_conf_cos_sim_matpltlib_sorted(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=None, rev=False)
    plot_velo_conf_cos_sim_matpltlib_sorted(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=None, rev=True)

    #plot_velo_conf_matpltlib_sorted(adata_total=total, dataset=dataset_short, method=method,fig_folder=fig_folder, split_seed=split_seed, celltype_label=None, rev=False)
    #plot_velo_conf_matpltlib_sorted(adata_total=total, dataset=dataset_short, method=method,fig_folder=fig_folder, split_seed=split_seed, celltype_label=None, rev=True)
    #plot_cosine_similarity_withRef_matpltlib_sorted(adata_split1=split1, adata_split2=split2, adata_total=total, dataset=dataset_short, method=method, fig_folder=fig_folder, split_seed=split_seed, rev=False)
    #plot_cosine_similarity_withRef_matpltlib_sorted(adata_split1=split1, adata_split2=split2, adata_total=total, dataset=dataset_short, method=method, fig_folder=fig_folder, split_seed=split_seed, rev=True)


plot_umap_color_sorted_pancreas(method='scv', split_seed=317)
plot_umap_color_sorted_pancreas(method='utv', split_seed=317)
plot_umap_color_sorted_pancreas(method='sct', split_seed=317)
plot_umap_color_sorted_pancreas(method='velovi', split_seed=317)
plot_umap_color_sorted_pancreas(method='velovi_woprep', split_seed=317)

plot_umap_color_sorted_pancreas(method='scv', split_seed=317, dataset_short='panINC', dataset_long = 'pancreasINC')
plot_umap_color_sorted_pancreas(method='utv', split_seed=317, dataset_short='panINC', dataset_long = 'pancreasINC')
plot_umap_color_sorted_pancreas(method='sct', split_seed=317, dataset_short='panINC', dataset_long = 'pancreasINC')
plot_umap_color_sorted_pancreas(method='velovi', split_seed=317, dataset_short='panINC', dataset_long = 'pancreasINC')
plot_umap_color_sorted_pancreas(method='velovi_woprep', split_seed=317, dataset_short='panINC', dataset_long = 'pancreasINC')

