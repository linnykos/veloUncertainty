import scvelo as scv
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib as plt

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

def get_dataset_long_and_short(dataset):
    if 'ery' in dataset:
        dataset_long = 'erythroid'
        dataset_short = 'ery'
    elif ('pan' in dataset) and ('INC' in dataset):
        dataset_long = 'pancreasINC'
        dataset_short = 'panINC'
    elif ('pan' in dataset) and (not 'INC' in dataset):
        dataset_long = 'pancreas'
        dataset_short = 'pan'
    elif ('la' in dataset) and ('M' in dataset):
        dataset_long = 'larryMult'
        dataset_short = 'larryMult'
    return dataset_long,dataset_short

def safe_len(data):
    if isinstance(data, np.ndarray):
        if data.ndim == 0:  return 0
        else: return data.size  # Return the number of elements in the array
    return 0  # For non-array types, return 0

def compute_cos_sim_rm_celltype_mean(s1,s2,dataset,method):
    celltype_label = get_celltype_label(dataset)
    res = []
    celltypes = s1.obs[celltype_label].cat.categories
    for ct in celltypes:
        ct_idx = np.where(s1.obs[celltype_label]==ct)[0]
        ct_s1 = s1.copy()[ct_idx,]
        ct_s2 = s2.copy()[ct_idx,]
        ct_s1.layers['velocity'] = ct_s1.layers['velocity']-np.mean(ct_s1.layers['velocity'], axis=0)
        ct_s2.layers['velocity'] = ct_s2.layers['velocity']-np.mean(ct_s2.layers['velocity'], axis=0)
        ct_cos_sim = compute_cosine_similarity_union(ct_s1,ct_s2,method)[0]
        ct_cos_sim = np.squeeze(ct_cos_sim)
        res.append(ct_cos_sim)
    return res


# global mean
def compute_cos_sim_rm_global_mean(s1,s2,method):
    s1_tmp = s1.copy()
    s2_tmp = s2.copy()
    vmean_s1 = np.mean(s1_tmp.layers['velocity'], axis=0)
    vmean_s2 = np.mean(s2_tmp.layers['velocity'], axis=0)
    s1_tmp.layers['velocity'] = s1_tmp.layers['velocity'] - vmean_s1
    s2_tmp.layers['velocity'] = s2_tmp.layers['velocity'] - vmean_s2
    return compute_cosine_similarity_union(s1_tmp,s2_tmp,method)

## plot
def plot_cos_sim_rm_celltype_mean_byCelltype(s1,s2,dataset,method,split_seed,fig_folder=None):
    dataset_long, dataset_short = get_dataset_long_and_short(dataset)
    if fig_folder is None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/'
        fig_folder += 'v4_'+dataset_long+'/seed'+str(split_seed)+'/cos_sim_rm_mean/'
    data_plot = compute_cos_sim_rm_celltype_mean(s1,s2,dataset_short,method)
    indices = [i for i, data in enumerate(data_plot) if safe_len(data) > 2]
    data_plot = [data_plot[i] for i in indices]
    celltype_label = get_celltype_label(dataset_short)
    celltypes = s1.obs[celltype_label].cat.categories[indices]
    plt.clf()
    plt.figure(figsize=(10, 6)) 
    plt.boxplot(data_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'),
               medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
               flierprops=dict(marker='o', color='rosybrown', alpha=0.5, markerfacecolor='aliceblue', markeredgecolor='rosybrown'))
    counts = [len(data) for data in data_plot]
    x_labels = [f'{celltypes[i]} (n={counts[i]})' for i in range(len(celltypes))]
    plt.xticks(ticks=range(1, len(celltypes) + 1), labels=x_labels, rotation=45, ha="right")
    plt.title(f'Cosine similarity with celltype-specific mean removed ({dataset_short} + {method})'+', split_seed='+str(split_seed))
    plt.ylim(-1,1)
    plt.tight_layout()
    plt.savefig(fig_folder+dataset_short+'_'+method+'_cos_sim_rm_celltype_mean_byCelltype_boxplot.png')
    # histograms
    plt.clf()
    n = len(data_plot)
    fig, axes = plt.subplots(int(np.ceil(n / 3)), 3, figsize=(12, 3 * np.ceil(n / 3)))
    axes = axes.ravel()
    for i, ct_data in enumerate(data_plot):
        axes[i].hist(ct_data, bins=30, alpha=0.7, color='lightsteelblue')
        axes[i].set_title(f'{celltypes[i]}'+', N='+str(len(ct_data)))
        axes[i].set_xlabel('Cosine Similarity')
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlim(-1,1)
    for j in range(i+1, len(axes)): fig.delaxes(axes[j])
    plt.suptitle('Cosine similarity with celltype-specific mean removed \n'+dataset_short+'+'+method+', seed='+str(split_seed))
    plt.tight_layout(rect=[0, 0, 1, 1]) 
    plt.savefig(fig_folder+dataset_short+'_'+method+'_cos_sim_rm_celltype_mean_byCelltype_hist.png')


def plot_cos_sim_rm_celltype_mean(s1,s2,dataset,method,split_seed,fig_folder=None):
    dataset_long, dataset_short = get_dataset_long_and_short(dataset)
    if fig_folder is None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/'
        fig_folder += 'v4_'+dataset_long+'/seed'+str(split_seed)+'/cos_sim_rm_mean/'
    data_plot = compute_cos_sim_rm_celltype_mean(s1,s2,dataset_short,method)
    indices = [i for i, data in enumerate(data_plot) if safe_len(data) > 2]
    data_plot = [data_plot[i] for i in indices]
    cos_sim = np.concatenate(data_plot).ravel()
    # histogram
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='gainsboro',color='powderblue') 
    max_frequency = np.max(counts)
    text_x = np.quantile(cos_sim,[.05])[0]
    text_y = max_frequency/5
    plt.axvline(np.mean(cos_sim), color='brown', linestyle='dashed', linewidth=1.5) ## add mean
    plt.axvline(np.median(cos_sim), color='peru', linestyle='dashed', linewidth=1.5) ## add median
    plt.text(text_x,text_y*2.5,'mean='+str(np.round(np.mean(cos_sim),4)), color='firebrick', fontsize=11)
    plt.text(text_x,text_y*3,'median='+str(np.round(np.median(cos_sim),4)), color='sienna', fontsize=11)
    plt.xlabel('cosine similarity')
    plt.ylabel('Frequency')
    plt.xlim(-1,1)
    plt.title('Cosine similarity with celltype-sepcific mean removed \n'+dataset_short+'+'+method+', split_seed='+str(split_seed))
    plt.savefig(fig_folder+dataset_short+'_'+method+'_cos_sim_rm_celltype_mean_allCells_hist.png')
    plt.clf()

def plot_cos_sim_rm_global_mean(s1,s2,dataset,method,split_seed,fig_folder=None):
    dataset_long, dataset_short = get_dataset_long_and_short(dataset)
    if fig_folder is None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/'
        fig_folder += 'v4_'+dataset_long+'/seed'+str(split_seed)+'/cos_sim_rm_mean/'
    cos_sim, Ngenes = compute_cos_sim_rm_global_mean(s1,s2,method)
    # histogram
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='gainsboro',color='powderblue') 
    max_frequency = np.max(counts)
    text_x = np.quantile(cos_sim,[.05])[0]
    text_y = max_frequency/5
    plt.axvline(np.mean(cos_sim), color='brown', linestyle='dashed', linewidth=1.5) ## add mean
    plt.axvline(np.median(cos_sim), color='peru', linestyle='dashed', linewidth=1.5) ## add median
    plt.text(text_x,text_y*2.5,'mean='+str(np.round(np.mean(cos_sim),4)), color='firebrick', fontsize=11)
    plt.text(text_x,text_y*3,'median='+str(np.round(np.median(cos_sim),4)), color='sienna', fontsize=11)
    plt.xlabel('cosine similarity')
    plt.ylabel('Frequency')
    plt.xlim(-1,1)
    plt.title('Cosine similarity with global mean removed, '+dataset+'+'+method+'\n Ngenes='+str(Ngenes)+', split_seed='+str(split_seed))
    plt.savefig(fig_folder+dataset+'_'+method+'_cos_sim_rm_global_mean_allCells_hist.png')
    plt.clf()

dataset_long = 'larryMult'
dataset_short = 'larryMult'
split_seed = 320
for method in ['scv','utv','sct','velovi','velovi_woprep']:
    outputAdded = True
    if method=='scv': outputAdded = False
    s1 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split1',allgenes=False,outputAdded=outputAdded)
    s2 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split2',allgenes=False,outputAdded=outputAdded)
    plot_cos_sim_rm_celltype_mean_byCelltype(s1,s2,dataset_short,method,split_seed,fig_folder=None)
    plot_cos_sim_rm_celltype_mean(s1,s2,dataset_short,method,split_seed,fig_folder=None)
    plot_cos_sim_rm_global_mean(s1,s2,dataset_short,method,split_seed,fig_folder=None)

def plot_cos_sim_rm_mean_allMethods(dataset,split_seed,methods=['scv','utv','sct','velovi','velovi_woprep']):
    dataset_long, dataset_short = get_dataset_long_and_short(dataset)
    for method in methods:
        outputAdded = True
        if method=='scv': outputAdded = False
        s1 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split1',allgenes=False,outputAdded=outputAdded)
        s2 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split2',allgenes=False,outputAdded=outputAdded)
        plot_cos_sim_rm_celltype_mean_byCelltype(s1,s2,dataset_short,method,split_seed,fig_folder=None)
        plot_cos_sim_rm_celltype_mean(s1,s2,dataset_short,method,split_seed,fig_folder=None)
        plot_cos_sim_rm_global_mean(s1,s2,dataset_short,method,split_seed,fig_folder=None)

plot_cos_sim_rm_mean_allMethods(dataset='panINC',split_seed=317)

