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



split_seed = 317

# erythroid
c_scv = get_cosine_similarity(method='scv', split_seed=split_seed, dataset='erythroid')
c_utv = get_cosine_similarity(method='utv', split_seed=split_seed, dataset='erythroid')
c_sct = get_cosine_similarity(method='sct', split_seed=split_seed, dataset='erythroid')
c_velovi = get_cosine_similarity(method='velovi', split_seed=split_seed, dataset='erythroid')
c_velov_w = get_cosine_similarity(method='velovi_woprep', split_seed=split_seed, dataset='erythroid')
data_to_plot = [c_scv, c_utv, c_sct, c_velovi, c_velov_w]
plot_cosine_similarity_boxplot_by_method(dataset='erythroid', split_seed=split_seed, data_to_plot=data_to_plot)

conf_scv = get_velocity_confidence(method='scv', dataset='erythroid')
conf_utv = get_velocity_confidence(method='utv', dataset='erythroid')
conf_sct = get_velocity_confidence(method='sct', dataset='erythroid')
conf_velovi = get_velocity_confidence(method='velovi', dataset='erythroid')
conf_velov_w = get_velocity_confidence(method='velovi_woprep', dataset='erythroid')
data_to_plot = [conf_scv, conf_utv, conf_sct, conf_velovi, conf_velov_w]
plot_velocity_confidence_boxplot_by_method(dataset='erythroid', data_to_plot=data_to_plot)


# pancreas
plot_cosine_similarity_boxplot_by_method(dataset='pancreas', split_seed=317, data_to_plot=None, fig_folder=None)
plot_velocity_confidence_boxplot_by_method(dataset='pancreas', split_seed=317, data_to_plot=None)

# pancreasINC
plot_cosine_similarity_boxplot_by_method(dataset='pancreasINC', split_seed=317, data_to_plot=None, fig_folder=None)
plot_velocity_confidence_boxplot_by_method(dataset='pancreasINC', split_seed=317, data_to_plot=None)

# larryMult
plot_velocity_confidence_boxplot_by_method(dataset='larryMult', split_seed=317, data_to_plot=None)

plot_cosine_similarity_boxplot_by_method(dataset='larryMult', split_seed=317, data_to_plot=None, fig_folder=None)
plot_cosine_similarity_boxplot_by_method(dataset='larryMult', split_seed=320, data_to_plot=None, fig_folder=None)
plot_cosine_similarity_boxplot_by_method(dataset='larryMult', split_seed=323, data_to_plot=None, fig_folder=None)
plot_cosine_similarity_boxplot_by_method(dataset='larryMult', split_seed=326, data_to_plot=None, fig_folder=None)
plot_cosine_similarity_boxplot_by_method(dataset='larryMult', split_seed=329, data_to_plot=None, fig_folder=None)


# plot gene correlation histogram

def get_dataset_names(dataset):
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
    return dataset_long, dataset_short
    
def plot_correlation_histogram(dataset, split_seed,fig_folder=None):
    dataset_long, dataset_short = get_dataset_names(dataset)
    if fig_folder==None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"
    split1_allgenes = read_data_v4(dataset_long,dataset_short,method=None,split_seed=split_seed,data_version='split1',allgenes=True,outputAdded=False)
    split2_allgenes = read_data_v4(dataset_long,dataset_short,method=None,split_seed=split_seed,data_version='split2',allgenes=True,outputAdded=False)
    split1_allgenes.layers['spliced_original'] = split1_allgenes.layers['spliced']
    split1_allgenes.layers['unspliced_original'] = split1_allgenes.layers['unspliced']
    split2_allgenes.layers['spliced_original'] = split2_allgenes.layers['spliced']
    split2_allgenes.layers['unspliced_original'] = split2_allgenes.layers['unspliced']
    common_genes = np.intersect1d(np.array(split1_allgenes.var.index), np.array(split2_allgenes.var.index))
    gene_names_split1 = split1_allgenes.var.index.copy()
    positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
    positions_split1 = [positions_dict_split1[gene] for gene in common_genes]
    gene_names_split2 = split2_allgenes.var.index.copy()
    positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
    positions_split2 = [positions_dict_split2[gene] for gene in common_genes]
    cor_spliced = compute_gene_correlation_between_splits(split1_allgenes.layers['spliced_original'][:,positions_split1],
                                                        split2_allgenes.layers['spliced_original'][:,positions_split2])
    cor_unspliced = compute_gene_correlation_between_splits(split1_allgenes.layers['unspliced_original'][:,positions_split1],
                                                            split2_allgenes.layers['unspliced_original'][:,positions_split2])
    Ngenes_spliced = len(cor_spliced[~np.isnan(cor_spliced)])
    Ngenes_unspliced = len(cor_unspliced[~np.isnan(cor_unspliced)])
    plt.clf()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
    axes[0].hist(cor_spliced, bins=20, edgecolor='gainsboro',color='powderblue') 
    axes[0].set_xlabel('Correlation')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Spliced splits gene correlation, '+dataset_long+'\n Ngenes='+str(Ngenes_spliced)+', split_seed='+str(split_seed))
    axes[1].hist(cor_unspliced, bins=20, edgecolor='gainsboro',color='powderblue') 
    axes[1].set_xlabel('Correlation')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Unspliced splits gene correlation, '+dataset_long+'\n Ngenes='+str(Ngenes_unspliced)+', split_seed='+str(split_seed))
    plt.tight_layout()
    plt.savefig(fig_folder+'gene_correlation_hist.png')
    plt.clf()


plot_correlation_histogram(dataset='ery', split_seed=317,fig_folder=None)
plot_correlation_histogram(dataset='pan', split_seed=317,fig_folder=None)
plot_correlation_histogram(dataset='panINC', split_seed=317,fig_folder=None)
plot_correlation_histogram(dataset='larryMult', split_seed=317,fig_folder=None)
plot_correlation_histogram(dataset='larryMult', split_seed=320,fig_folder=None)
plot_correlation_histogram(dataset='larryMult', split_seed=323,fig_folder=None)
plot_correlation_histogram(dataset='larryMult', split_seed=326,fig_folder=None)
plot_correlation_histogram(dataset='larryMult', split_seed=329,fig_folder=None)


def plot_overdispersion_histogram(dataset, split_seed,fig_folder=None):
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
    dataset_long, dataset_short = get_dataset_names(dataset)
    if fig_folder==None:
        fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"
    overdisp_S = np.array(pd.read_csv(data_folder+'v4_'+dataset_long+'/'+dataset_short+'_overdisp_S.csv')['x'])
    overdisp_U = np.array(pd.read_csv(data_folder+'v4_'+dataset_long+'/'+dataset_short+'_overdisp_U.csv')['x'])
    S_xlab = 'Overdispersion'
    U_xlab = 'Overdispersion'
    if np.max(overdisp_S)>100: S_xlab = S_xlab+' (capped by 100)'
    if np.max(overdisp_U)>100: U_xlab = U_xlab+' (capped by 100)'
    plt.clf()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
    axes[0].hist(overdisp_S.clip(max=100), bins=30, edgecolor='gainsboro',color='powderblue') 
    axes[0].set_xlabel(S_xlab)
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Spliced splits gene overdispersion, '+dataset_long+'\n Ngenes='+str(np.sum(~np.isnan(overdisp_S)))+', split_seed='+str(split_seed))
    axes[1].hist(overdisp_U.clip(max=100), bins=30, edgecolor='gainsboro',color='powderblue') 
    axes[1].set_xlabel(U_xlab)
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Unspliced splits gene overdispersion, '+dataset_long+'\n Ngenes='+str(np.sum(~np.isnan(overdisp_U)))+', split_seed='+str(split_seed))
    plt.tight_layout()
    plt.savefig(fig_folder+'overdispersion_hist.png')
    plt.clf()

plot_overdispersion_histogram(dataset='ery', split_seed=317,fig_folder=None)
plot_overdispersion_histogram(dataset='pan', split_seed=317,fig_folder=None)
plot_overdispersion_histogram(dataset='panINC', split_seed=317,fig_folder=None)
plot_overdispersion_histogram(dataset='larryMult', split_seed=317,fig_folder=None)

