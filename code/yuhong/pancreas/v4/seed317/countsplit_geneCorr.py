import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import scvelo as scv

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

split_seed = 317
dataset_long = 'pancreas'
dataset_short = 'pan'
celltype_label = 'clusters'


data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/' 

# correlation between splits
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

total = read_raw_adata(dataset_short)
total.layers['spliced_original'] = total.layers['spliced'].copy()
total.layers['unspliced_original'] = total.layers['unspliced'].copy()
celltypes = np.array(list(total.obs[celltype_label].cat.categories))

df = pd.DataFrame()
df['cor_spliced'] = cor_spliced
df['cor_unspliced'] = cor_unspliced


##### corr_fracNonzero_allgenes
n_rows = total.layers['spliced_original'].shape[0]
nonzeros_per_column_S = total.layers['spliced_original'].getnnz(axis=0)
fraction_nonzeros_S = nonzeros_per_column_S / n_rows
df['frac_nnz_S'] = fraction_nonzeros_S
nonzeros_per_column_U = total.layers['unspliced_original'].getnnz(axis=0)
fraction_nonzeros_U = nonzeros_per_column_U / n_rows
df['frac_nnz_U'] = fraction_nonzeros_U

plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(df['frac_nnz_S'], df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title('Correlation between splits (seed='+str(split_seed)+', spliced), N='+str(len(df['cor_spliced'][~np.isnan(df['cor_spliced'])])))
axes[0].set_xlabel('fraction of nonzeros (Spliced)')
axes[0].set_ylabel('Correlation')
axes[1].scatter(df['frac_nnz_U'], df['cor_unspliced'], color='royalblue',alpha=0.4)
axes[1].set_title('Correlation between splits (seed='+str(split_seed)+',unspliced), N='+str(len(df['cor_unspliced'][~np.isnan(df['cor_unspliced'])])))
axes[1].set_xlabel('fraction of nonzeros (Unspliced)')
axes[1].set_ylabel('Correlation')
plt.tight_layout()
plt.savefig(fig_folder+'corr_fracNonzero_allgenes.png') 

