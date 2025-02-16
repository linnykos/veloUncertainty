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

indicator_GPC = genes.isin(genes_GPC)

## we would like to control for the density of each gene across cells
#adata_GPC = adata[:, genes_GPC]
#adata_nGPC = adata[:, genes_nGPC]

## The likelihood computation requires recover_dynamics and is not available in stochastic or deterministic mode.
## to compute gene R2, it is unavoidable to filter out a few genes
fnzS_all = np.asarray(np.mean(adata.layers["spliced"] > 0, axis=0))[0]
fnzU_all = np.asarray(np.mean(adata.layers["unspliced"] > 0, axis=0))[0]

scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode='deterministic')
r2_all = np.array(adata.var.velocity_r2)


df = pd.DataFrame( {'GPC': indicator_GPC, 'fnzS': fnzS_all, 'fnzU': fnzU_all, 'r2': r2_all} )
df.to_csv(data_folder+'v4_greenleaf/glf_r2_df.csv')


"""
fnzS_GPC = fnzS_all[indicator_GPC]
fnzU_GPC = fnzU_all[indicator_GPC]
r2_GPC = r2_all[indicator_GPC]

fnzS_nGPC = fnzS_all[~indicator_GPC]
fnzU_nGPC = fnzU_all[~indicator_GPC]
r2_nGPC = r2_all[~indicator_GPC]

df_GPC = pd.DataFrame({'fnzS': fnzS_GPC, 'fnzU': fnzU_GPC, 'r2': r2_GPC })
df_nGPC = pd.DataFrame({'fnzS': fnzS_nGPC, 'fnzU': fnzU_nGPC, 'r2': r2_nGPC })
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors
import pandas as pd

dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

df = pd.read_csv(data_folder+'v4_greenleaf/glf_r2_df.csv')
df_GPC = df[df['GPC']]
df_nGPC = df[~df['GPC']]

# Parameters
bs = 4  # Adjust as needed
w = 1.0  # Standard deviation for normal distribution

# Step 1: Generate bins
bins1 = np.linspace(np.min(df_GPC['fnzS']), np.max(df_GPC['fnzS']), bs)
bins2 = np.linspace(np.min(df_GPC['fnzU']), np.max(df_GPC['fnzU']), bs)
bins3 = np.linspace(np.min(df_GPC['r2']), np.max(df_GPC['r2']), bs)
bin_data = np.array([[b1, b2, b3] for b1 in bins1 for b2 in bins2 for b3 in bins3]) # Create bin_data (similar to R's `do.call(rbind, ...)`)
# Step 2: Compute Euclidean distances
bin_dist = cdist(bin_data, bin_data, metric='euclidean')  # Pairwise distances
# Step 3: Compute bin probabilities using normal distribution
bin_p = norm.pdf(bin_dist, loc=0, scale=w) # dim: (bs^3)x(bs^3)
# Step 4: Compute bin membership using kNN (k=1)
nbrs = NearestNeighbors(n_neighbors=1).fit(bin_data)
membership_GPC = nbrs.kneighbors(df_GPC[['fnzS','fnzU','r2']], return_distance=False).flatten()
membership_nGPC = nbrs.kneighbors(df_nGPC[['fnzS','fnzU','r2']], return_distance=False).flatten()

# the frequency table of GPC genes' membership
unique, counts = np.unique(membership_GPC, return_counts=True)
freq_GPC = dict(zip(unique, counts))

n_sample_scale = 2
genes_select = np.array([])
np.random.seed(227)
for lab in freq_GPC.keys():
    print(lab)
    n_sample = freq_GPC[lab]*n_sample_scale
    idx_nGPC = (np.where( membership_nGPC==lab )[0])
    idx_nGPC = idx_nGPC[~np.isin(idx_nGPC, genes_select)]
    if (len(idx_nGPC) <= n_sample):
        genes_select = np.append( genes_select, idx_nGPC ) 
        n_sample -= len(idx_nGPC)
        while (n_sample>0):
            print('n_sample='+str(n_sample))
            val = bin_dist[lab,]
            val[lab] = np.Inf
            lab = (np.random.choice( np.where(val==np.min(val))[0], size=1 )[0]) # sample a nbr bin to be sampled from, if there are multiple
            idx_nGPC = (np.where( membership_nGPC==lab )[0])
            idx_nGPC = idx_nGPC[~np.isin(idx_nGPC, genes_select)]
            if len(idx_nGPC)==0:
                print('->' + str(lab)+', case 1')
                continue
            elif len(idx_nGPC) <= n_sample:
                print('->' + str(lab)+', case 2')
                genes_select = np.append(genes_select, idx_nGPC)
                n_sample -= len(idx_nGPC)
            else:
                print('->' + str(lab)+', case 3')
                idx_sample = np.random.choice(idx_nGPC, size=n_sample, replace=False, p=None)
                genes_select = np.append(genes_select, idx_sample)
                n_sample = 0
    else:
        idx_sample = np.random.choice(idx_nGPC, size=n_sample, replace=False, p=None)
        genes_select = np.append(genes_select, idx_sample)
    #print('N(genes=0): '+str(np.sum(genes_select==0))+', Nduplicated: '+str(len(genes_select) - len(np.unique(genes_select))))

genes_select = np.unique(genes_select)   

import scanpy as sc
adata = sc.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 
genes = adata.var.index
genes_GPC_all = pd.read_csv(data_folder+'v4_'+dataset_long+'/glf_GPC.csv')
genes_GPC = np.intersect1d(adata.var.index, genes_GPC_all) # 184 genes
genes_nGPC = genes[~genes.isin(genes_GPC)] # 32464 genes

adata_nGPC = adata[:, genes_nGPC]
adata_nGPC_subset = adata_nGPC[:, [int(x) for x in genes_select]]

split_seed = 317
split1 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
split2 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

split1_nGPC = split1[:, genes_nGPC]
split2_nGPC = split2[:, genes_nGPC]
split1_nGPC_subset = split1_nGPC[:, [int(x) for x in genes_select]]
split2_nGPC_subset = split2_nGPC[:, [int(x) for x in genes_select]]

adata_nGPC_subset.var['highly_variable'] = True
split1_nGPC_subset.var['highly_variable'] = True
split2_nGPC_subset.var['highly_variable'] = True
# cd /home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/
adata_nGPC_subset.write(data_folder+'v4_greenleaf/glf_total_nGPCgrid.h5ad')
split1_nGPC_subset.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split1_nGPCgrid.h5ad')
split2_nGPC_subset.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split2_nGPCgrid.h5ad')

"""
# Step 5: Compute bin density (equivalent to tabulate2)
density_GPC = np.bincount(membership_GPC, minlength=bs**3)
density_nGPC = np.bincount(membership_nGPC, minlength=bs**3)
"""



