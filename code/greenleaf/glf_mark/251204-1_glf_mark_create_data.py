import pandas as pd

dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'
genes_mark = pd.read_csv(data_folder+'marker_genes_top25.csv')
genes_mark = genes_mark['names'].unique()


import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

genes_GPC = pd.read_csv(data_folder+'glf_GPC.csv')



#genes_gsea_df = pd.read_csv('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/glf_genes_gsea/glf_genes_gsea.csv')
#genes_gsea = genes_gsea_df['ENSEMBL']

total = ad.read_h5ad(data_folder+'glf_total_allgenes.h5ad') 

np.intersect1d(genes_mark, genes_GPC) # 45 genes
genes_mark_selected = np.intersect1d(genes_mark, total.var.index) # 


# seed317
split_seed = 317
split1 = sc.read(data_folder+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
split2 = sc.read(data_folder+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

total[:, genes_mark_selected].write(data_folder+'/glf_genes_mark/glf_total_mark_genes.h5ad')
#split1[:, genes_mark_selected].write(data_folder+'/glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split1_mark_genes.h5ad')
#split2[:, genes_mark_selected].write(data_folder+'/glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split2_mark_genes.h5ad')

# more seeds
for split_seed in [317,320,323,326,329]:
    split1 = sc.read(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
    split2 = sc.read(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')
    split1[:, genes_mark_selected].write(data_folder+'glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split1_mark_genes.h5ad')
    split2[:, genes_mark_selected].write(data_folder+'glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split2_mark_genes.h5ad')
    print(str(split_seed)+' done!')



def select_genes_blk(freq_GPC, n_sample_scale, grid_seed=227):
    np.random.seed(grid_seed)
    genes_select = np.array([])
    for lab in freq_GPC.keys():
        print(lab)
        n_sample = freq_GPC[lab]*n_sample_scale
        idx_nGPC = (np.where( membership_ctrl==lab )[0])
        idx_nGPC = idx_nGPC[~np.isin(idx_nGPC, genes_select)]
        if (len(idx_nGPC) <= n_sample):
            genes_select = np.append( genes_select, idx_nGPC ) 
            n_sample -= len(idx_nGPC)
            nbr_list = np.array([lab]) 
            while (n_sample>0):
                print('n_sample='+str(n_sample))
                val = bin_dist[lab,]
                val[nbr_list] = np.Inf
                lab_nbr = (np.random.choice( np.where(val==np.min(val))[0], size=1 )[0])
                nbr_list = np.append(nbr_list, lab_nbr)
                idx_nGPC = (np.where( membership_ctrl==lab_nbr )[0])
                idx_nGPC = idx_nGPC[~np.isin(idx_nGPC, genes_select)]
                if len(idx_nGPC)==0:
                    print('->' + str(lab_nbr)+', case 1')
                    continue
                elif len(idx_nGPC) <= n_sample:
                    print('->' + str(lab_nbr)+', case 2')
                    genes_select = np.append(genes_select, idx_nGPC)
                    n_sample -= len(idx_nGPC)
                else:
                    print('->' + str(lab_nbr)+', case 3, len(idx_nGPC)>n_sample:', str(len(idx_nGPC)>n_sample))
                    idx_sample = np.random.choice(idx_nGPC, size=n_sample, replace=False, p=None)
                    genes_select = np.append(genes_select, idx_sample)
                    n_sample = 0
        else:
            idx_sample = np.random.choice(idx_nGPC, size=n_sample, replace=False, p=None)
            genes_select = np.append(genes_select, idx_sample)
        #print('N(genes=0): '+str(np.sum(genes_select==0))+', Nduplicated: '+str(len(genes_select) - len(np.unique(genes_select))))
    return genes_select

def get_adata_nGPCblk(genes_nGPC, freq_GPC, adata_nGPC, split_seed, grid_seed, n_sample_scale):
    genes_select = select_genes_blk(freq_GPC, n_sample_scale=n_sample_scale, grid_seed=grid_seed)
    adata_nGPC_subset = adata_nGPC[:, [int(x) for x in genes_select]]
    split1 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
    split2 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')
    split1_nGPC = split1[:, genes_nGPC]
    split2_nGPC = split2[:, genes_nGPC]
    split1_nGPC_subset = split1_nGPC[:, [int(x) for x in genes_select]]
    split2_nGPC_subset = split2_nGPC[:, [int(x) for x in genes_select]]
    adata_nGPC_subset.var['highly_variable'] = True
    split1_nGPC_subset.var['highly_variable'] = True
    split2_nGPC_subset.var['highly_variable'] = True
    return adata_nGPC_subset, split1_nGPC_subset, split2_nGPC_subset


import requests
def get_ensembl_id(gene_name):
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    response = requests.get(url)
    if response.ok:
        data = response.json()
        return data[0]["id"] if data else None
    return None


# 0. setup
dataset_long = 'greenleaf'
dataset_short = 'glf'
#data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
genes_mark = pd.read_csv(data_folder+'marker_genes_top25.csv')
genes_mark = genes_mark['names'].unique()

adata = ad.read_h5ad(data_folder+'glf_total_allgenes.h5ad') 
genes_mark_selected = np.intersect1d(genes_mark, adata.var.index) # 

#adata = sc.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 
genes = adata.var.index
#genes_GPC_all = pd.read_csv(data_folder+'v4_'+dataset_long+'/glf_GPC.csv')
#genes_GPC = np.intersect1d(adata.var.index, genes_GPC_all) # 184 genes

"""
genes_blk = pd.read_csv(data_folder+'v4_greenleaf/Writeup13_blacklist_mark_gpc_genes.csv', header=None)
ensembl_ids = {gene: get_ensembl_id(gene) for gene in genes_blk.values.flatten()}.values()
pd.DataFrame({'gene_id': ensembl_ids, 'gene_name': genes_blk.values.flatten()}).to_csv(data_folder+'v4_greenleaf/genes_blacklist_mark_gpc.csv')
genes_blk_ids = np.array(list( {gene: get_ensembl_id(gene) for gene in genes_blk.values.flatten()}.values() ), dtype=genes_GPC.dtype)
"""

genes_GPC = pd.read_csv(data_folder+'glf_GPC.csv')

genes_blk_ids = pd.read_csv(data_folder+'genes_blacklist_gsea_gpc.csv')['gene_id']
genes_ctrl = genes[~genes.isin( np.union1d(genes_GPC, np.union1d(genes_blk_ids, genes_mark) ) )] # 
adata_ctrl = adata[:, genes_ctrl]

## R2
indicator_mark = genes.isin( genes_mark )
indicator_ctrl = genes.isin( genes_ctrl )
## The likelihood computation requires recover_dynamics and is not available in stochastic or deterministic mode.
## to compute gene R2, it is unavoidable to filter out a few genes
fnzS_all = np.asarray(np.mean(adata.layers["spliced"] > 0, axis=0))[0]
fnzU_all = np.asarray(np.mean(adata.layers["unspliced"] > 0, axis=0))[0]

import scvelo as scv
adata_r2 = adata.copy()
scv.pp.normalize_per_cell(adata_r2)
scv.pp.log1p(adata_r2)
scv.pp.moments(adata_r2, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_r2, mode='deterministic')
r2_all = np.array(adata_r2.var.velocity_r2)

df = pd.DataFrame( {'gene_id':adata_r2.var.index, 'ismark': indicator_mark, 'ctrl': indicator_ctrl, 
					'fnzS': fnzS_all, 'fnzU': fnzU_all, 'r2': r2_all} )
df.to_csv(data_folder+'glf_r2_df_mark.csv')

df = pd.read_csv(data_folder+'glf_r2_df_mark.csv')
df_mark = df[df['ismark']]
df_ctrl = df[df['ctrl']]

from scipy.spatial.distance import cdist
from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors

# Parameters
bs = 4  # Adjust as needed
w = 1.0  # Standard deviation for normal distribution
# Step 1: Generate bins
bins1 = np.linspace(np.min(df_mark['fnzS']), np.max(df_mark['fnzS']), bs)
bins2 = np.linspace(np.min(df_mark['fnzU']), np.max(df_mark['fnzU']), bs)
bins3 = np.linspace(np.min(df_mark['r2']), np.max(df_mark['r2']), bs)
bin_data = np.array([[b1, b2, b3] for b1 in bins1 for b2 in bins2 for b3 in bins3]) 
# Step 2: Compute Euclidean distances
bin_dist = cdist(bin_data, bin_data, metric='euclidean')  # Pairwise distances
# Step 3: Compute bin probabilities using normal distribution
bin_p = norm.pdf(bin_dist, loc=0, scale=w) # dim: (bs^3)x(bs^3)
# Step 4: Compute bin membership using kNN (k=1)
nbrs = NearestNeighbors(n_neighbors=1).fit(bin_data)
membership_mark = nbrs.kneighbors(df_mark[['fnzS','fnzU','r2']], return_distance=False).flatten()
membership_ctrl = nbrs.kneighbors(df_ctrl[['fnzS','fnzU','r2']], return_distance=False).flatten()
# the frequency table of GPC genes' membership
unique, counts = np.unique(membership_mark, return_counts=True)
freq_mark = dict(zip(unique, counts))

gene_set_name = 'markctrl'
n_sample_scale = 2

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
for i in range(5):
    split_seed = [317,320,323,326,329][i]
    grid_seed = [227,230,233,236,239][i]
    total_nGPC, split1_nGPC, split2_nGPC = get_adata_nGPCblk(genes_nGPC=genes_ctrl, freq_GPC=freq_mark, adata_nGPC=adata_ctrl, 
                                                             split_seed=split_seed, grid_seed=grid_seed, n_sample_scale=n_sample_scale)
    total_nGPC.write(data_folder+'v4_greenleaf/glf_genes_mark/glf_total_'+gene_set_name+str(grid_seed)+'.h5ad')
    split1_nGPC.write(data_folder+'v4_greenleaf/glf_genes_mark/glf_split1_'+gene_set_name+str(grid_seed)+'.h5ad')
    split2_nGPC.write(data_folder+'v4_greenleaf/glf_genes_mark/glf_split2_'+gene_set_name+str(grid_seed)+'.h5ad')
    print('##################### '+str(split_seed)+' done')





