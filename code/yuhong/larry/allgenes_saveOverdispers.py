import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import scvelo as scv
from scipy.stats import gamma, beta
from scipy.sparse import csr_matrix, find
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import compute_gene_correlation_between_splits

def estimate_overdisps(X):
    res = []
    p = None
    if X.ndim==1: 
        p = 1
        data = { 'counts': X }
        df = pd.DataFrame(data)
        model = smf.negativebinomial('counts ~ 1', data=df)
        result = model.fit()
        res.append(1/result.params['alpha'])
    if X.ndim==2:
        p = X.shape[1]
        for col in range(p):
            print(col)
            y = X[:, col]
            if np.sum(y)==0:
                res.append(np.inf)
            else: 
                if hasattr(y, 'todense'):
                    y = y.todense().A 
                df = pd.DataFrame({'counts': y.flatten()})
                model = smf.negativebinomial('counts ~ 1', data=df)
                result = model.fit()
                alpha = result.params['alpha']
                b = 1/alpha
                res.append(b)
    return np.array(res)


data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

print('************************* Read data *************************')
total = ad.read_h5ad(data_folder+"v2_larry/larry.h5ad")
split1_allgenes = sc.read_h5ad(data_folder+'v2_larry/larry_split1_allgenes.h5ad')
split2_allgenes = sc.read_h5ad(data_folder+'v2_larry/larry_split2_allgenes.h5ad')

print('************************* Read data done *************************')
celltype_label = 'state_info'
celltypes = np.array(list(total.obs[celltype_label].cat.categories))
total_S = total.layers['spliced']
total_U = total.layers['unspliced']

print('************************* Calculating correlations *************************')
cor_spliced = compute_gene_correlation_between_splits(split1_allgenes.layers['spliced'],
                                                      split2_allgenes.layers['spliced'])
cor_unspliced = compute_gene_correlation_between_splits(split1_allgenes.layers['unspliced'],
                                                        split2_allgenes.layers['unspliced'])

print('************************* Estimating overdispersion parameters *************************')
overdisps_S = estimate_overdisps(total_S)
overdisps_U = estimate_overdisps(total_U)

print('************************* Estimating DE *************************')
scv.pp.normalize_per_cell(total)
scv.pp.log1p(total)
sc.tl.rank_genes_groups(total, groupby=celltype_label, method="wilcoxon")
df = pd.DataFrame(columns=celltypes)
for ct in celltypes:
    pval = sc.get.rank_genes_groups_df(total, group=ct)['pvals']
    df[ct] = pval

df['cor_spliced'] = cor_spliced
df['cor_unspliced'] = cor_unspliced
df['overdisps_S'] = overdisps_S
df['overdisps_U'] = overdisps_U

df.to_csv(data_folder+'v2_larry/larry_df_allgenes.csv')

df['gene_names'] = total.var.index

def compute_nnz(adata,type):
    n_rows = adata.layers[type].shape[0]
    nonzeros_per_column = adata.layers[type].getnnz(axis=0)
    fraction_nonzeros = nonzeros_per_column / n_rows
    return fraction_nonzeros

df['frac_nnz_S'] = compute_nnz(total,'spliced')
df['frac_nnz_U'] = compute_nnz(total,'unspliced')


print('************************* Writing csv *************************')
df.to_csv(data_folder+'v2_larry/larry_df_allgenes.csv')
print('************************* All done *************************')