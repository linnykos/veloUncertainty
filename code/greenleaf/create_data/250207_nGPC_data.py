import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
import anndata as ad

dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
adata = ad.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 

genes = adata.var.index
genes_GPC_all = pd.read_csv(data_folder+'v4_'+dataset_long+'/glf_GPC.csv')
genes_GPC = np.intersect1d(adata.var.index, genes_GPC_all) # 184 genes
genes_nGPC = genes[~genes.isin(genes_GPC)] # 32464 genes

## we would like to control for the density of each gene across cells
# adata[:, adata.var.index.isin(genes_GPC)]
adata_GPC = adata[:, genes_GPC]
adata_nGPC = adata[:, genes_nGPC]

### sparsity of all genes
np.array(np.mean(adata.layers['spliced']>0, 0 ))[0]

# fraction of nonzero in spliced/unspliced counts, for GPC and non-GPC genes
fnzS_GPC = np.asarray(np.mean(adata_GPC.layers["spliced"] > 0, axis=0)).reshape(-1, 1)
fnzU_GPC = np.asarray(np.mean(adata_GPC.layers["unspliced"] > 0, axis=0)).reshape(-1, 1)

fnzS_nGPC = np.asarray(np.mean(adata_nGPC.layers["spliced"] > 0, axis=0)).reshape(-1, 1)
fnzU_nGPC = np.asarray(np.mean(adata_nGPC.layers["unspliced"] > 0, axis=0)).reshape(-1, 1)

# data frame containing the sparsity/nonzero fraction information
fnz_GPC = pd.DataFrame({'S': np.asarray(np.mean(adata_GPC.layers["spliced"] > 0, axis=0))[0], 
                        'U': np.asarray(np.mean(adata_GPC.layers["unspliced"] > 0, axis=0))[0]})
fnz_nGPC = pd.DataFrame({'S': np.asarray(np.mean(adata_nGPC.layers["spliced"] > 0, axis=0))[0], 
                         'U': np.asarray(np.mean(adata_nGPC.layers["unspliced"] > 0, axis=0))[0]})

########################################################
## k-means
def get_freq_table(kmeans_lab, prop=False):
    unique, counts = np.unique(kmeans_lab, return_counts=True)
    print( dict(zip(unique, counts)) )
    if (prop):
        counts = np.round( counts/len(kmeans_lab), 4 )
        print( dict(zip(unique, counts)) )
        return dict(zip(unique, counts))
    else: return dict(zip(unique, counts))

from sklearn.cluster import KMeans
kmeans_GPC = KMeans(n_clusters=10, random_state=0, n_init="auto").fit( fnz_GPC )
kmeans_GPC.labels_
freq_GPC = get_freq_table(kmeans_GPC.labels_)

# use the GPC genes' clusters to predict non-GPC genes
kmeans_nGPC = kmeans_GPC.predict(fnz_nGPC)
get_freq_table(kmeans_nGPC)

np.random.seed(227)
n_sample_scale = 2
genes_select = np.array([])
for lab in freq_GPC.keys():
    n_sample = freq_GPC[lab]*n_sample_scale
    idx_nGPC = np.where(kmeans_nGPC==lab)[0]
    idx_sample = np.random.choice(idx_nGPC, size=n_sample, replace=False, p=None)
    genes_select = np.append(genes_select, idx_sample)

adata_nGPC_subset = adata_nGPC[:, [int(x) for x in genes_select]]


split_seed = 317
split1 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
split2 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

split1_nGPC = split1[:, genes_nGPC]
split2_nGPC = split2[:, genes_nGPC]
split1_nGPC_subset = split1_nGPC[:, [int(x) for x in genes_select]]
split2_nGPC_subset = split2_nGPC[:, [int(x) for x in genes_select]]

# cd /home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/
adata_nGPC_subset.write(data_folder+'v4_greenleaf/glf_total_nGPC227.h5ad')
split1_nGPC_subset.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_nGPC227.h5ad')
split2_nGPC_subset.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_nGPC227.h5ad')

"""
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *


import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'scv_nGPC227'
split_seed=317
scv_compute_velocity(split1_nGPC_subset,dataset_short) 

plot_velocity_scv_utv(adata_in=split1_nGPC_subset,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
"""