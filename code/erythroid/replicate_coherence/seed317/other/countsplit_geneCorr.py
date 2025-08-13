import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import scvelo as scv

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import compute_gene_correlation_between_splits

split_seed = 317
dataset_long = 'erythroid'
dataset_short = 'ery'
celltype_label = 'celltype'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
#fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/' 

# correlation between splits
split1_allgenes = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
split2_allgenes = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

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

total = sc.read_h5ad(data_folder+'Gastrulation/erythroid_lineage.h5ad')
total.layers['spliced_original'] = total.layers['spliced'].copy()
total.layers['unspliced_original'] = total.layers['unspliced'].copy()
celltypes = np.array(list(total.obs[celltype_label].cat.categories))


scv.pp.normalize_per_cell(total)
scv.pp.log1p(total)
sc.tl.rank_genes_groups(total, groupby=celltype_label, method='wilcoxon')
pval_df = pd.DataFrame(columns=celltypes)
for ct in celltypes:
    pval = sc.get.rank_genes_groups_df(total, group=ct)['pvals']
    pval_df[ct] = pval
    
pval_df['min'] = pval_df.min(axis=1) # must be computed right after the pvalue loop
pval_df['cor_spliced'] = cor_spliced
pval_df['cor_unspliced'] = cor_unspliced

pval_df['-log10min'] = -np.log10(pval_df['min'])
pval_df['-log10min'] = pval_df['-log10min'].clip(upper=30)


##### corr_pvals_allgenes
plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_df['-log10min'], pval_df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title('Gene correlation between splits ('+'seed='+str(split_seed)+',spliced), N='+str(Ngenes_spliced))
axes[0].set_xlabel('-log10(pval)')
axes[0].set_ylabel('Correlation')
# Plotting unspliced
axes[1].scatter(pval_df['-log10min'], pval_df['cor_unspliced'],color='royalblue',alpha=0.4)
axes[1].set_title('Gene correlation between splits ('+'seed='+str(split_seed)+',unspliced), N='+str(Ngenes_unspliced))
axes[1].set_xlabel('-log10(pval)')
axes[1].set_ylabel('Correlation')
# Adjusting layout to avoid overlap
plt.tight_layout()
plt.savefig(fig_folder+'corr_pvals_allgenes1.png') 

##### corr_overdisps_allgenes
overdisp_S = np.array(pd.read_csv(data_folder+'v4_'+dataset_long+'/'+dataset_short+'_overdisp_S.csv')['x'])
overdisp_U = np.array(pd.read_csv(data_folder+'v4_'+dataset_long+'/'+dataset_short+'_overdisp_U.csv')['x'])
# these are in the order of genes in total

(np.sum(overdisp_S == np.inf),np.sum(overdisp_S != np.inf)) # (42090, 11711)
(np.sum(overdisp_U == np.inf),np.sum(overdisp_U != np.inf)) # (44225, 9576)

overdisp_S[overdisp_S == np.inf] = np.nan
overdisp_U[overdisp_U == np.inf] = np.nan

np.quantile(overdisp_S[~np.isnan(overdisp_S)],[0.,.25,.5,.75,1.])
# array([8.10151656e-06, 4.49853838e-01, 2.11927380e+00, 4.75772047e+00, 2.41311903e+06])

pval_df['overdisps_S'] = overdisp_S
pval_df['overdisps_S'] = pval_df['overdisps_S'].clip(upper=50)

pval_df['overdisps_U'] = overdisp_U
pval_df['overdisps_U'] = pval_df['overdisps_U'].clip(upper=50)

#np.sum(~np.isnan(pval_df['overdisps_S']+pval_df['cor_spliced']))

plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
axes[0].scatter(pval_df['overdisps_S'], pval_df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title('Correlation between splits (seed=317, spliced), N='+str(np.sum(~np.isnan(pval_df['overdisps_S']+pval_df['cor_spliced']))))
axes[0].set_xlabel('overdispersion')
axes[0].set_ylabel('Correlation')
axes[1].scatter(pval_df['overdisps_U'], pval_df['cor_unspliced'],color='royalblue',alpha=0.4)
axes[1].set_title('Correlation between splits (seed=317, unspliced), N='+str(np.sum(~np.isnan(pval_df['overdisps_U']+pval_df['cor_unspliced']))))
axes[1].set_xlabel('overdispersion')
axes[1].set_ylabel('Correlation')
# Adjusting layout to avoid overlap
plt.tight_layout()
plt.savefig(fig_folder+'corr_overdisps_allgenes.png') 

## clip overdispersion
plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_df['overdisps_S'].clip(upper=100), pval_df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title('Correlation between splits (seed=317, spliced), N='+str(np.sum(~np.isnan(pval_df['overdisps_S']+pval_df['cor_spliced']))))
axes[0].set_xlabel('overdispersion')
axes[0].set_ylabel('Correlation')
axes[1].scatter(pval_df['overdisps_U'].clip(upper=100), pval_df['cor_unspliced'],color='royalblue',alpha=0.4)
axes[1].set_title('Correlation between splits (seed=317, unspliced), N='+str(np.sum(~np.isnan(pval_df['overdisps_U']+pval_df['cor_unspliced']))))
axes[1].set_xlabel('overdispersion')
axes[1].set_ylabel('Correlation')
# Adjusting layout to avoid overlap
plt.tight_layout()
plt.savefig(fig_folder+'corr_overdisps_clip100_allgenes.png') 



##### corr_fracNonzero_allgenes
#total2 = ad.read_h5ad(data_folder+'v4_larry/larry.h5ad')
n_rows = total.layers['spliced_original'].shape[0]
nonzeros_per_column_S = total.layers['spliced_original'].getnnz(axis=0)
fraction_nonzeros_S = nonzeros_per_column_S / n_rows
pval_df['frac_nnz_S'] = fraction_nonzeros_S
nonzeros_per_column_U = total.layers['unspliced_original'].getnnz(axis=0)
fraction_nonzeros_U = nonzeros_per_column_U / n_rows
pval_df['frac_nnz_U'] = fraction_nonzeros_U

plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_df['frac_nnz_S'], pval_df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title('Correlation between splits (seed='+str(split_seed)+', spliced), N='+str(len(pval_df['cor_spliced'][~np.isnan(pval_df['cor_spliced'])])))
axes[0].set_xlabel('fraction of nonzeros (Spliced)')
axes[0].set_ylabel('Correlation')
axes[1].scatter(pval_df['frac_nnz_U'], pval_df['cor_unspliced'], color='royalblue',alpha=0.4)
axes[1].set_title('Correlation between splits (seed='+str(split_seed)+',unspliced), N='+str(len(pval_df['cor_unspliced'][~np.isnan(pval_df['cor_unspliced'])])))
axes[1].set_xlabel('fraction of nonzeros (Unspliced)')
axes[1].set_ylabel('Correlation')
plt.tight_layout()
plt.savefig(fig_folder+'corr_fracNonzero_allgenes.png') 


########################################
import matplotlib.pyplot as plt

df = pval_df

pval_min = np.array((df[['Blood progenitors 1','Blood progenitors 2','Erythroid1','Erythroid2','Erythroid3']]).min(axis=1))
pval_mean = np.array((df[['Blood progenitors 1','Blood progenitors 2','Erythroid1','Erythroid2','Erythroid3']]).mean(axis=1))
pval_median = np.array((df[['Blood progenitors 1','Blood progenitors 2','Erythroid1','Erythroid2','Erythroid3']]).median(axis=1))

pval_log_min = np.array((-np.log10(df[['Blood progenitors 1','Blood progenitors 2','Erythroid1','Erythroid2','Erythroid3']])).clip(upper=60).min(axis=1))
pval_log_mean = np.array((-np.log10(df[['Blood progenitors 1','Blood progenitors 2','Erythroid1','Erythroid2','Erythroid3']])).clip(upper=60).mean(axis=1))
pval_log_median = np.array((-np.log10(df[['Blood progenitors 1','Blood progenitors 2','Erythroid1','Erythroid2','Erythroid3']])).clip(upper=60).median(axis=1))

# RuntimeWarning: divide by zero encountered in log10

Ngenes_spliced = np.sum(~np.isnan(df['cor_spliced']))
Ngenes_unspliced = np.sum(~np.isnan(df['cor_unspliced']))


plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_log_mean, df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title('Correlation of gene expr between splits (spliced), N='+str(Ngenes_spliced))
axes[0].set_xlabel('mean -log10(pval)')
axes[0].set_ylabel('Correlation')
# Plotting unspliced
axes[1].scatter(pval_log_mean, df['cor_unspliced'],color='royalblue',alpha=0.4)
axes[1].set_title('Correlation of gene expr between splits (unspliced), N='+str(Ngenes_unspliced))
axes[1].set_xlabel('mean -log10(pval)')
axes[1].set_ylabel('Correlation')
plt.tight_layout()
plt.savefig(fig_folder+'corr_pvals_mean_neglog10_allgenes.png') 


plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(pval_log_median, df['cor_spliced'], color='royalblue',alpha=0.4)
axes[0].set_title('Correlation of gene expr between splits (spliced), N='+str(Ngenes_spliced))
axes[0].set_xlabel('median -log10(pval)')
axes[0].set_ylabel('Correlation')
# Plotting unspliced
axes[1].scatter(pval_log_median, df['cor_unspliced'],color='royalblue',alpha=0.4)
axes[1].set_title('Correlation of gene expr between splits (unspliced), N='+str(Ngenes_unspliced))
axes[1].set_xlabel('median -log10(pval)')
axes[1].set_ylabel('Correlation')
# Adjusting layout to avoid overlap
plt.tight_layout()
plt.savefig(fig_folder+'corr_pvals_median_neglog10_allgenes.png') 


## overdispersion ~ nonzero%
plt.clf()
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
# Plotting spliced
axes[0].scatter(df['frac_nnz_S'][~np.isnan(df['cor_spliced'])], df['overdisps_S'][~np.isnan(df['cor_spliced'])], color='royalblue',alpha=0.4)
axes[0].set_title('Overdispersion vs nonzero% (seed=317,spliced), N='+str(np.sum(~np.isnan(df['cor_spliced']))))
axes[0].set_xlabel('Fraction of nonzeros (Spliced)')
axes[0].set_ylabel('Overdispersion parameter')
axes[1].scatter(df['frac_nnz_U'][~np.isnan(df['cor_unspliced'])], df['overdisps_U'][~np.isnan(df['cor_unspliced'])], color='royalblue',alpha=0.4)
axes[1].set_title('Overdispersion vs nonzero% (seed=317,unspliced), N='+str(np.sum(~np.isnan(df['cor_unspliced']))))
axes[1].set_xlabel('Fraction of nonzeros (Unspliced)')
axes[1].set_ylabel('Overdispersion parameter')
plt.savefig(fig_folder+'overdispersion_fracNonzero_allgenes.png') 

