import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import datetime

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/" 

total = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_total_v2.h5ad')
split1 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split1_v2.h5ad')
split2 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split2_v2.h5ad')

### common genes being filtered out
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index))
common_genes.shape # 1481 common genes
np.intersect1d(np.array(split1.var.index), np.array(total.var.index)).shape # 1335
np.intersect1d(np.array(split2.var.index), np.array(total.var.index)).shape # 1328

# get the counts of commonly selected genes
gene_names_split1 = split1.var.index.copy()
positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
positions_split1 = [positions_dict_split1[gene] for gene in common_genes]

gene_names_split2 = split2.var.index.copy()
positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
positions_split2 = [positions_dict_split2[gene] for gene in common_genes]

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

# split1.layers['spliced'][:,positions_split1] is cell-by-gene, (9815, 1237)
cor = compute_gene_correlation_between_splits(split1.layers['spliced_original'][:,positions_split1],split2.layers['spliced_original'][:,positions_split2])
cor[~np.isnan(cor)].shape
np.quantile(cor[~np.isnan(cor)],[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

## plot correlation between splits (all in one color)
def plot_gene_correlation_between_splits(adata1,adata2,fig_path,split_seed=317):
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
    # Plotting the first dataset
    axes[0].scatter(range(Ngenes_spliced), cor_spliced[~np.isnan(cor_spliced)],color='royalblue',alpha=0.4)
    axes[0].set_title("Correlation of gene expr between splits (spliced), N="+str(Ngenes_spliced))
    axes[0].set_xlabel("Genes")
    axes[0].set_ylabel("Correlation")
    # Plotting the second dataset
    axes[1].scatter(range(Ngenes_unspliced), cor_unspliced[~np.isnan(cor_unspliced)],color='royalblue',alpha=0.4)
    axes[1].set_title("Correlation of gene expr between splits (unspliced), N="+str(Ngenes_unspliced))
    axes[1].set_xlabel("Genes")
    axes[1].set_ylabel("Correlation")
    # Adjusting layout to avoid overlap
    plt.tight_layout()
    plt.savefig(fig_folder+fig_path) 
    plt.clf()

plot_gene_correlation_between_splits(split1,split2,fig_path="scv/split_cor/corr_between_splits_colorSame.png")

# color by sparsity

# color by housekeeping genes
