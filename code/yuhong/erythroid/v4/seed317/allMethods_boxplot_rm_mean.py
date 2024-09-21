import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# celltype mean
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

def get_cosine_similarity_rm_celltype_mean(method, split_seed, dataset):
    dataset_long,dataset_short = get_dataset_long_and_short(dataset)
    outputAdded=True
    if method=='scv' or method=='utv': outputAdded=False
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=outputAdded)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=outputAdded)
    cos_sim = compute_cos_sim_rm_celltype_mean(split1,split2,dataset,method)
    indices = [i for i, data in enumerate(cos_sim) if safe_len(data) > 2]
    cos_sim = [cos_sim[i] for i in indices]
    cos_sim = np.concatenate(cos_sim).ravel()
    return cos_sim

# global mean
def compute_cos_sim_rm_global_mean(s1,s2,method):
    s1_tmp = s1.copy()
    s2_tmp = s2.copy()
    vmean_s1 = np.mean(s1_tmp.layers['velocity'], axis=0)
    vmean_s2 = np.mean(s2_tmp.layers['velocity'], axis=0)
    s1_tmp.layers['velocity'] = s1_tmp.layers['velocity'] - vmean_s1
    s2_tmp.layers['velocity'] = s2_tmp.layers['velocity'] - vmean_s2
    return compute_cosine_similarity_union(s1_tmp,s2_tmp,method)

def get_cosine_similarity_rm_global_mean(method, split_seed, dataset):
    dataset_long,dataset_short = get_dataset_long_and_short(dataset)
    outputAdded=True
    if method=='scv' or method=='utv': outputAdded=False
    split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=outputAdded)
    split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=outputAdded)
    cos_sim, Ngenes = compute_cos_sim_rm_global_mean(split1,split2,method)
    return cos_sim

def get_cosine_similarity_rm_mean(mean_type, method, split_seed, dataset):
    if 'cell' in mean_type: return get_cosine_similarity_rm_celltype_mean(method, split_seed, dataset)
    elif 'glo' in mean_type: return get_cosine_similarity_rm_global_mean(method, split_seed, dataset)

                                  
def plot_cosine_similarity_rm_mean_boxplot_by_method(mean_type, dataset, split_seed, data_to_plot=None, fig_folder=None):
    dataset_long,dataset_short = get_dataset_long_and_short(dataset)
    if fig_folder==None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"
    if data_to_plot==None:
        c_scv = get_cosine_similarity_rm_mean(mean_type=mean_type,method='scv', split_seed=split_seed, dataset=dataset)
        c_utv = get_cosine_similarity_rm_mean(mean_type=mean_type,method='utv', split_seed=split_seed, dataset=dataset)
        c_sct = get_cosine_similarity_rm_mean(mean_type=mean_type,method='sct', split_seed=split_seed, dataset=dataset)
        c_velovi = get_cosine_similarity_rm_mean(mean_type=mean_type,method='velovi', split_seed=split_seed, dataset=dataset)
        c_velov_w = get_cosine_similarity_rm_mean(mean_type=mean_type,method='velovi_woprep', split_seed=split_seed, dataset=dataset)
        data_to_plot = [c_scv, c_utv, c_sct, c_velovi, c_velov_w]
    plt.clf()
    plt.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'),
                medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
                flierprops=dict(marker='o', color='rosybrown', alpha=0.5, markerfacecolor='aliceblue', markeredgecolor='rosybrown'))
    plt.title('Velocity cosine similarity for '+dataset+'across methods with '+mean_type+' mean removed, seed='+str(split_seed))
    plt.xticks([1, 2, 3, 4, 5], ['scv', 'utv', 'sct', 'velovi', 'velovi_woprep'])
    plt.ylim(-1, 1)
    plt.savefig(fig_folder+'cos_sim_rm_'+mean_type+'_mean_byMethod_boxplot.png')
    plt.clf()




# pancreas
plot_cosine_similarity_rm_mean_boxplot_by_method(mean_type='global',dataset='pancreas', split_seed=317, data_to_plot=None, fig_folder=None)
plot_cosine_similarity_rm_mean_boxplot_by_method(mean_type='celltype',dataset='pancreas', split_seed=317, data_to_plot=None, fig_folder=None)

# erythroid
plot_cosine_similarity_rm_mean_boxplot_by_method(mean_type='global',dataset='ery', split_seed=317, data_to_plot=None, fig_folder=None)
plot_cosine_similarity_rm_mean_boxplot_by_method(mean_type='celltype',dataset='ery', split_seed=317, data_to_plot=None, fig_folder=None)


