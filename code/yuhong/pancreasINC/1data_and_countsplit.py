import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

split_seed = 317

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

adata = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
cell_index = np.array(np.where(adata.obs['clusters']!='Pre-endocrine')[0])

"""
adata.var # 'highly_variable_genes'
adata.uns # 'clusters_coarse_colors', 'clusters_colors', 'day_colors', 'neighbors', 'pca'
adata.obsm # 'X_pca', 'X_umap'
adata.layers # 'spliced', 'unspliced'
adata.obsp # 'distances', 'connectivities'

# keep adata.var and uns unchanged
adata.obs = adata.obs[adata.obs['clusters']!='Pre-endocrine']
adata.obsm['X_pca'] = adata.obsm['X_pca'][cell_index,]
adata.obsm['X_umap'] = adata.obsm['X_umap'][cell_index,]
adata.layers['spliced'] = adata.layers['spliced'][cell_index,:]
adata.layers['unspliced'] = adata.layers['unspliced'][cell_index,:]
del adata.obsp
"""

S_mat = adata.layers['spliced'][cell_index,:].copy()
U_mat = adata.layers['unspliced'][cell_index,:].copy()

def create_adata_pancreasINC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}#adata_old.uns['clusters_colors'].copy()
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    del adata_new.obs['S_score']
    del adata_new.obs['G2M_score']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

total = create_adata_pancreasINC(S=S_mat,U=U_mat,adata_old=adata)
total.write(data_folder+'v2_pancreasINC/pancreasINC_total_allgenes.h5ad')

def run_countsplit_with_overdispersion(S,U,split_seed):
    overdisps_S = estimate_overdisps(S)
    overdisps_U = estimate_overdisps(U)
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisps_S)
    u1, u2  = countsplit(U,overdisps=overdisps_U)
    return [[s1,u1],[s2,u2]]#[s1,s2,u1,u2]

def create_adata_pancreas(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.uns = {'clusters_colors':adata.uns['clusters_colors'].copy()}
    adata_split.obs = pd.DataFrame(index=adata_total.obs.index)
    for obs_col in adata_total.obs.columns:
        adata_split.obs[obs_col] = adata_total.obs[obs_col].copy()
    adata_split.var = pd.DataFrame(index=adata_total.var.index)
    for var_col in adata_total.var.columns:
        adata_split.var[var_col] = adata_total.var[var_col].copy()
    return adata_split

def countsplit_and_create_adata(S,U,total,split_seed):
    split1,split2 = run_countsplit_with_overdispersion(S=S,U=U,split_seed=split_seed)
    counts_adata1 = ad.AnnData(X=split1[0].astype(np.float32))
    counts_adata1.layers["spliced"] = split1[0]
    counts_adata1.layers["unspliced"] = split1[1]
    counts_adata2 = ad.AnnData(X=split2[0].astype(np.float32))
    counts_adata2.layers["spliced"] = split2[0]
    counts_adata2.layers["unspliced"] = split2[1]
    adata1 = create_adata_pancreas(split1[0],split1[1],total)
    adata2 = create_adata_pancreas(split2[0],split2[1],total)
    return adata1,adata2

adata_split1,adata_split2 = countsplit_and_create_adata(S=S_mat,U=U_mat,total=total,split_seed=317)

adata_split1.write(data_folder+'v2_pancreasINC/pancreasINC_split1_allgenes.h5ad')
adata_split2.write(data_folder+'v2_pancreasINC/pancreasINC_split2_allgenes.h5ad')

from v2_functions import compute_gene_correlation_between_splits,plot_gene_correlation_between_splits
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

cor_spliced = compute_gene_correlation_between_splits(adata_split1.layers['spliced'],adata_split2.layers['spliced'])
"""
Number of valid values = 15734
Quantiles: [-1.64251512e-01 -1.58638966e-02 -1.02339088e-02 -6.22471493e-03
 -3.23854538e-03 -1.44285724e-03 -5.58365860e-04  2.67208308e-03
  9.36738426e-03  1.93597664e-02  9.96119127e-01]
"""
cor_unspliced = compute_gene_correlation_between_splits(adata_split1.layers['unspliced'],adata_split2.layers['unspliced'])
"""
Number of valid values = 14231
Quantiles: [-1.45782025e-01 -1.39385092e-02 -8.83841825e-03 -5.60111177e-03
 -3.28255081e-03 -1.70776287e-03 -7.89775842e-04 -3.22268772e-04
  7.02317667e-03  1.99176679e-02  9.10826485e-01]
"""

fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_pancreasINC/'

## plot correlation between splits (all in one color)
plt.clf()
plt.scatter(range(len(cor_spliced[~np.isnan(cor_spliced)])), cor_spliced[~np.isnan(cor_spliced)],color='royalblue',alpha=0.5)
plt.title("Correlation of gene expr between splits (spliced counts, all genes)")
plt.xlabel("genes (with nonzero expr)")
plt.ylabel("correlation")
plt.savefig(fig_folder+"corr_between_splits_allgenes_spliced.png") 
plt.clf()

plt.clf()
plt.scatter(range(len(cor_unspliced[~np.isnan(cor_unspliced)])), cor_unspliced[~np.isnan(cor_unspliced)],color='royalblue',alpha=0.5)
plt.title("Correlation of gene expr between splits (spliced counts, all genes)")
plt.xlabel("genes (with nonzero expr)")
plt.ylabel("correlation")
plt.savefig(fig_folder+"corr_between_splits_allgenes_unspliced.png") 
plt.clf()
