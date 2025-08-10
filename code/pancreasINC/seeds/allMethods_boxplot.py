import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

def get_velocity_confidence(method, dataset):
    if 'ery' in dataset:
        dataset_long = 'erythroid'
        dataset_short = 'ery'
    elif 'pan' in dataset and (not 'INC' in dataset):
        dataset_long = 'pancreas'
        dataset_short = 'pan'
    elif 'pan' in dataset and ('INC' in dataset):
        dataset_long = 'pancreasINC'
        dataset_short = 'panINC'
    elif 'larry' in dataset and ('M' in dataset):
        dataset_long = 'larryMult'
        dataset_short = 'larryMult'
    if not method=='scv':
        total = read_data_v4(dataset_long,dataset_short,method,split_seed=317,data_version='total',allgenes=False,outputAdded=True)
    else:
        total = read_data_v4(dataset_long,dataset_short,method,split_seed=317,data_version='total',allgenes=False,outputAdded=False)
    if not 'velocity_confidence' in total.obs.columns:
        scv.tl.velocity_confidence(total)
    conf = total.obs['velocity_confidence']
    return conf

def plot_velocity_confidence_boxplot_by_method(dataset, split_seed, data_to_plot=None, fig_folder=None):
    if fig_folder==None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset+'/seed'+str(split_seed)+"/"
    if data_to_plot==None:
        c_scv = get_velocity_confidence(method='scv',dataset=dataset)
        c_utv = get_velocity_confidence(method='utv',dataset=dataset)
        c_sct = get_velocity_confidence(method='sct',dataset=dataset)
        c_velovi = get_velocity_confidence(method='velovi',dataset=dataset)
        c_velov_w = get_velocity_confidence(method='velovi_woprep',dataset=dataset)
        data_to_plot = [c_scv, c_utv, c_sct, c_velovi, c_velov_w]
    plt.clf()
    plt.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'),
                medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
                flierprops=dict(marker='o', color='rosybrown', alpha=0.5, markerfacecolor='aliceblue', markeredgecolor='rosybrown'))
    plt.title('Velocity confidence for '+dataset+' across methods')
    plt.xticks([1, 2, 3, 4, 5], ['scv', 'utv', 'sct', 'velovi', 'velovi_woprep'])
    plt.ylim(-1, 1)
    plt.savefig(fig_folder+'conf_byMethod_boxplot.png')
    plt.clf()


def get_cosine_similarity(method, split_seed, dataset):
    if 'ery' in dataset:
        dataset_long = 'erythroid'
        dataset_short = 'ery'
    elif 'pan' in dataset and (not 'INC' in dataset):
        dataset_long = 'pancreas'
        dataset_short = 'pan'
    elif 'pan' in dataset and ('INC' in dataset):
        dataset_long = 'pancreasINC'
        dataset_short = 'panINC'
    elif 'larry' in dataset and ('M' in dataset):
        dataset_long = 'larryMult'
        dataset_short = 'larryMult'
    if method=='scv' or method=='utv':
        split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
        split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
    else:
        split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
        split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
    cos_sim = compute_cosine_similarity_union(split1,split2,method)[0]
    return cos_sim

def plot_cosine_similarity_boxplot_by_method(dataset, split_seed, data_to_plot=None, fig_folder=None):
    if fig_folder==None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset+'/seed'+str(split_seed)+"/"
    if data_to_plot==None:
        c_scv = get_cosine_similarity(method='scv', split_seed=split_seed, dataset=dataset)
        c_utv = get_cosine_similarity(method='utv', split_seed=split_seed, dataset=dataset)
        c_sct = get_cosine_similarity(method='sct', split_seed=split_seed, dataset=dataset)
        c_velovi = get_cosine_similarity(method='velovi', split_seed=split_seed, dataset=dataset)
        c_velov_w = get_cosine_similarity(method='velovi_woprep', split_seed=split_seed, dataset=dataset)
        data_to_plot = [c_scv, c_utv, c_sct, c_velovi, c_velov_w]
    plt.clf()
    plt.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'),
                medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
                flierprops=dict(marker='o', color='rosybrown', alpha=0.5, markerfacecolor='aliceblue', markeredgecolor='rosybrown'))
    plt.title('Velocity cosine similarity for '+dataset+' across methods, seed='+str(split_seed))
    plt.xticks([1, 2, 3, 4, 5], ['scv', 'utv', 'sct', 'velovi', 'velovi_woprep'])
    plt.ylim(-1, 1)
    plt.savefig(fig_folder+'cos_sim_byMethod_boxplot.png')
    plt.clf()



split_seed = 320
plot_cosine_similarity_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None, fig_folder=None)
plot_velocity_confidence_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None)

split_seed = 323
plot_cosine_similarity_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None, fig_folder=None)
plot_velocity_confidence_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None)

split_seed = 326
plot_cosine_similarity_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None, fig_folder=None)
plot_velocity_confidence_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None)

split_seed = 329
plot_cosine_similarity_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None, fig_folder=None)
plot_velocity_confidence_boxplot_by_method(dataset='pancreasINC', split_seed=split_seed, data_to_plot=None)
