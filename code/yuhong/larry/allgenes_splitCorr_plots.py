import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import scvelo as scv

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import compute_gene_correlation_between_splits


data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_larry/" 

# correlation between splits
split1_allgenes = sc.read_h5ad(data_folder+'v2_larry/larry_split1_allgenes.h5ad')
split2_allgenes = sc.read_h5ad(data_folder+'v2_larry/larry_split2_allgenes.h5ad')

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

total = ad.read_h5ad(data_folder+"v2_larry/larry.h5ad")
celltype_label = 'state_info'
celltypes = np.array(list(total.obs[celltype_label].cat.categories))

scv.pp.normalize_per_cell(total)
scv.pp.log1p(total)
sc.tl.rank_genes_groups(total, groupby=celltype_label, method="wilcoxon")
pval_df = pd.DataFrame(columns=celltypes)
for ct in celltypes:
    pval = sc.get.rank_genes_groups_df(total, group=ct)['pvals'][positions_split1]
    pval_df[ct] = pval
    
pval_df['min'] = pval_df.min(axis=1) # must be computed right after the pvalue loop
pval_df['cor_spliced'] = cor_spliced
pval_df['cor_unspliced'] = cor_unspliced

pval_df['-log10min'] = -np.log10(pval_df['min'])
pval_df['-log10min'] = pval_df['-log10min'].clip(upper=30)

np.max(1,-np.log10(pval_df.min(axis=1).round(5)))

##### corr_pvals_allgenes
plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_df['-log10min'], pval_df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title("Correlation of gene expr between splits (spliced), N="+str(Ngenes_spliced))
axes[0].set_xlabel("-log10(pval)")
axes[0].set_ylabel("Correlation")
# Plotting unspliced
axes[1].scatter(pval_df['-log10min'], pval_df['cor_unspliced'],color='royalblue',alpha=0.4)
axes[1].set_title("Correlation of gene expr between splits (unspliced), N="+str(Ngenes_unspliced))
axes[1].set_xlabel("-log10(pval)")
axes[1].set_ylabel("Correlation")
# Adjusting layout to avoid overlap
plt.tight_layout()
plt.savefig(fig_folder+'corr_pvals_allgenes.png') 

##### corr_overdisps_allgenes
from countsplit import estimate_overdisps
total2 = ad.read_h5ad(data_folder+"v2_larry/larry.h5ad")
overdisps_S = estimate_overdisps(total2.layers['spliced'])
overdisps_U = estimate_overdisps(total2.layers['unspliced'])

overdisps_S[overdisps_S == np.inf] = np.nan
overdisps_U[overdisps_U == np.inf] = np.nan

pval_df['overdisps_S'] = overdisps_S[positions_split1]
pval_df['overdisps_S'] = pval_df['overdisps_S'].clip(upper=50)

pval_df['overdisps_U'] = overdisps_U[positions_split1]
pval_df['overdisps_U'] = pval_df['overdisps_U'].clip(upper=50)

#np.sum(~np.isnan(pval_df['overdisps_S']+pval_df['cor_spliced']))

plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_df['overdisps_S'], pval_df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title("Correlation of gene expr between splits (spliced), N="+str(np.sum(~np.isnan(pval_df['overdisps_S']+pval_df['cor_spliced']))))
axes[0].set_xlabel("overdispersion")
axes[0].set_ylabel("Correlation")
axes[1].scatter(pval_df['overdisps_U'], pval_df['cor_unspliced'],color='royalblue',alpha=0.4)
axes[1].set_title("Correlation of gene expr between splits (unspliced), N="+str(np.sum(~np.isnan(pval_df['overdisps_U']+pval_df['cor_unspliced']))))
axes[1].set_xlabel("overdispersion")
axes[1].set_ylabel("Correlation")
# Adjusting layout to avoid overlap
plt.tight_layout()
plt.savefig(fig_folder+'corr_overdisps_allgenes.png') 


##### corr_fracNonzero_allgenes
n_rows = total2.layers['spliced'].shape[0]
nonzeros_per_column = total2.layers['spliced'].getnnz(axis=0)
fraction_nonzeros = nonzeros_per_column / n_rows
pval_df['frac_nnz'] = fraction_nonzeros[positions_split1]

plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_df['frac_nnz'], pval_df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title("Correlation of gene expr between splits (spliced), N="+str(len(pval_df['cor_spliced'][~np.isnan(pval_df['cor_spliced'])])))
axes[0].set_xlabel("fraction of nonzeros")
axes[0].set_ylabel("Correlation")
axes[1].scatter(pval_df['frac_nnz'], pval_df['cor_unspliced'], color='royalblue',alpha=0.4)
axes[1].set_title("Correlation of gene expr between splits (spliced), N="+str(len(pval_df['cor_unspliced'][~np.isnan(pval_df['cor_unspliced'])])))
axes[1].set_xlabel("fraction of nonzeros")
axes[1].set_ylabel("Correlation")
plt.tight_layout()
plt.savefig(fig_folder+'corr_fracNonzero_allgenes.png') 



########################################
## sct
dataset_long = 'larry'
dataset_short = 'larry'
method = 'sct'

total=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_split2_v2.h5ad')

common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index))
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
Ngenes_unspliced = len(cor_unspliced[~np.isnan(cor_unspliced)])

scv.pp.normalize_per_cell(total)
scv.pp.log1p(total)
sc.tl.rank_genes_groups(total, groupby=celltype_label, method="wilcoxon")
celltype_pval_min = {}
celltype_pval_median = {}
for ct in celltypes:
    pval = sc.get.rank_genes_groups_df(total, group=ct)['pvals']
    celltype_pval_min[ct] = np.min(pval)
    celltype_pval_median[ct] = np.median(pval)


len(np.intersect1d(np.array(split1.var.index), np.array(total.var.index)))
len(np.intersect1d(np.array(total.var.index), np.array(split2.var.index)))
len(np.intersect1d(np.array(total.var.index), common_genes))
# different sets of genes in total, split1, and split2
