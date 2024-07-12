import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

method = 'velovi'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/" 

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 9815 × 310
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 9815 × 119
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 9815 × 128

### common genes being filtered out
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) # 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 91
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(split1.var.index),np.array(total.var.index)).shape[0])) # 114
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(split2.var.index),np.array(total.var.index)).shape[0])) # 125

def compute_gene_correlation_between_splits(mat1,mat2):
    if hasattr(mat1, 'todense'):
        mat1 = mat1.todense().A 
        mat2 = mat2.todense().A 
    m1 = np.transpose(mat1)
    m2 = np.transpose(mat2)
    cor = []
    for i in range(m1.shape[0]):
        if i%1000==0:
            print(i)
        if np.sum(m1[i,:])==0 and np.sum(m2[i,:])==0:
            cor.append(np.nan)
        else:
            cor.append(np.corrcoef(m1[i,:],m2[i,:])[0,1])
    cor = np.array(cor)
    print("Number of valid values = "+str(cor[~np.isnan(cor)].shape[0]))
    print("Quantiles: "+str(np.quantile(cor[~np.isnan(cor)],[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])))
    return cor

def plot_gene_correlation_between_splits(adata1,adata2,fig_path):
    common_genes = np.intersect1d(np.array(adata1.var.index), np.array(adata2.var.index))
    gene_names_split1 = split1.var.index.copy()
    positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
    positions_split1 = [positions_dict_split1[gene] for gene in common_genes]
    gene_names_split2 = split2.var.index.copy()
    positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
    positions_split2 = [positions_dict_split2[gene] for gene in common_genes]
    cor_spliced = compute_gene_correlation_between_splits(split1.layers['spliced_original'][:,positions_split1],
                                                          split2.layers['spliced_original'][:,positions_split2])
    cor_unspliced = compute_gene_correlation_between_splits(split1.layers['unspliced_original'][:,positions_split1],
                                                            split2.layers['unspliced_original'][:,positions_split2])
    Ngenes_spliced = len(cor_spliced[~np.isnan(cor_spliced)])
    Ngenes_unspliced = len(cor_spliced[~np.isnan(cor_spliced)])
    plt.clf()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
    # Plotting spliced
    axes[0].scatter(range(Ngenes_spliced), cor_spliced[~np.isnan(cor_spliced)],color='royalblue',alpha=0.4)
    axes[0].set_title("Correlation of gene expr between splits (spliced), N="+str(Ngenes_spliced))
    axes[0].set_xlabel("Genes")
    axes[0].set_ylabel("Correlation")
    # Plotting unspliced
    axes[1].scatter(range(Ngenes_unspliced), cor_unspliced[~np.isnan(cor_unspliced)],color='royalblue',alpha=0.4)
    axes[1].set_title("Correlation of gene expr between splits (unspliced), N="+str(Ngenes_unspliced))
    axes[1].set_xlabel("Genes")
    axes[1].set_ylabel("Correlation")
    # Adjusting layout to avoid overlap
    plt.tight_layout()
    plt.savefig(fig_folder+fig_path) 
    plt.clf()

plot_gene_correlation_between_splits(split1,split2,fig_path=method+"/split_cor/corr_between_splits_colorSame.png")

# color by sparsity

# color by housekeeping genes




### common genes being selected and used for velocity computation
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index))
common_genes.shape #  common genes
np.intersect1d(np.array(split1.var.index), np.array(total.var.index)).shape # 
np.intersect1d(np.array(split2.var.index), np.array(total.var.index)).shape # 

# get the counts of commonly selected genes
gene_names_split1 = split1.var.index.copy()
positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
positions_split1 = [positions_dict_split1[gene] for gene in common_genes]

gene_names_split2 = split2.var.index.copy()
positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
positions_split2 = [positions_dict_split2[gene] for gene in common_genes]
