import scanpy as sc
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors
import pandas as pd


## functions
def select_genes(n_sample_scale, grid_seed=227):
    np.random.seed(grid_seed)
    genes_select = np.array([])
    for lab in freq_GPC.keys():
        print(lab)
        n_sample = freq_GPC[lab]*n_sample_scale
        idx_nGPC = (np.where( membership_nGPC==lab )[0])
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
                idx_nGPC = (np.where( membership_nGPC==lab_nbr )[0])
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

def get_adata_nGPC(adata_nGPC, split_seed, grid_seed, n_sample_scale):
    genes_select = select_genes(n_sample_scale=n_sample_scale, grid_seed=grid_seed)
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

# 0. setup
dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

adata = sc.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 
genes = adata.var.index
genes_GPC_all = pd.read_csv(data_folder+'v4_'+dataset_long+'/glf_GPC.csv')
genes_GPC = np.intersect1d(adata.var.index, genes_GPC_all) # 184 genes
genes_nGPC = genes[~genes.isin(genes_GPC)] # 32464 genes
adata_nGPC = adata[:, genes_nGPC]

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
bin_data = np.array([[b1, b2, b3] for b1 in bins1 for b2 in bins2 for b3 in bins3]) 
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


# 1. select genes using countsplit_seed=317, grid_seed=227
split_seed = 317
grid_seed = 227
n_sample_scale = 2
total_nGPC, split1_nGPC, split2_nGPC = get_adata_nGPC(adata_nGPC=adata_nGPC, split_seed=split_seed, grid_seed=grid_seed, n_sample_scale=n_sample_scale)

total_nGPC.write(data_folder+'v4_greenleaf/glf_total_nGPCgrid'+str(grid_seed)+'.h5ad')
split1_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split1_nGPCgrid'+str(grid_seed)+'.h5ad')
split2_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split2_nGPCgrid'+str(grid_seed)+'.h5ad')

genes227 = total_nGPC.var.indices

# 2. select genes using countsplit_seed=320, grid_seed=230
split_seed = 320
grid_seed = 230
n_sample_scale = 2
total_nGPC, split1_nGPC, split2_nGPC = get_adata_nGPC(adata_nGPC=adata_nGPC, split_seed=split_seed, grid_seed=grid_seed, n_sample_scale=n_sample_scale)

total_nGPC.write(data_folder+'v4_greenleaf/glf_total_nGPCgrid'+str(grid_seed)+'.h5ad')
split1_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split1_nGPCgrid'+str(grid_seed)+'.h5ad')
split2_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split2_nGPCgrid'+str(grid_seed)+'.h5ad')

genes230 = total_nGPC.var.indices

# 3. select genes using countsplit_seed=323, grid_seed=233
split_seed = 323
grid_seed = 233
n_sample_scale = 2
total_nGPC, split1_nGPC, split2_nGPC = get_adata_nGPC(adata_nGPC=adata_nGPC, split_seed=split_seed, grid_seed=grid_seed, n_sample_scale=n_sample_scale)

total_nGPC.write(data_folder+'v4_greenleaf/glf_total_nGPCgrid'+str(grid_seed)+'.h5ad')
split1_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split1_nGPCgrid'+str(grid_seed)+'.h5ad')
split2_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split2_nGPCgrid'+str(grid_seed)+'.h5ad')

genes233= total_nGPC.var.indices

# 4. select genes using countsplit_seed=326, grid_seed=236
split_seed = 326
grid_seed = 236
n_sample_scale = 2
total_nGPC, split1_nGPC, split2_nGPC = get_adata_nGPC(adata_nGPC=adata_nGPC, split_seed=split_seed, grid_seed=grid_seed, n_sample_scale=n_sample_scale)

total_nGPC.write(data_folder+'v4_greenleaf/glf_total_nGPCgrid'+str(grid_seed)+'.h5ad')
split1_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split1_nGPCgrid'+str(grid_seed)+'.h5ad')
split2_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split2_nGPCgrid'+str(grid_seed)+'.h5ad')

genes236 = total_nGPC.var.indices

# 5. select genes using countsplit_seed=329, grid_seed=239
split_seed = 329
grid_seed = 239
n_sample_scale = 2
total_nGPC, split1_nGPC, split2_nGPC = get_adata_nGPC(adata_nGPC=adata_nGPC, split_seed=split_seed, grid_seed=grid_seed, n_sample_scale=n_sample_scale)

total_nGPC.write(data_folder+'v4_greenleaf/glf_total_nGPCgrid'+str(grid_seed)+'.h5ad')
split1_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split1_nGPCgrid'+str(grid_seed)+'.h5ad')
split2_nGPC.write(data_folder+'v4_'+dataset_long+'/glf_seed'+str(split_seed)+'_split2_nGPCgrid'+str(grid_seed)+'.h5ad')

genes239 = total_nGPC.var.indices

# fix the problem that total.X are normalized
adata = sc.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 
for i in range(5):
    split_seed = [317,320,323,326,329][i]
    grid_seed = [227,230,233,236,239][i]
    genes_total = sc.read(data_folder+'v4_greenleaf/glf_total_nGPCgrid'+str(grid_seed)+'.h5ad').var.index
    adata[:,genes_total].write(data_folder+'v4_greenleaf/glf_total_nGPCgrid'+str(grid_seed)+'.h5ad')


# save the control genes
df_genes_ctrl = pd.DataFrame({'genes227': np.array(sc.read(data_folder+'v4_'+dataset_long+'/glf_seed317_split1_nGPCgrid227.h5ad').var.index),
              'genes230': np.array(sc.read(data_folder+'v4_'+dataset_long+'/glf_seed320_split1_nGPCgrid230.h5ad').var.index ),
              'genes233': np.array(sc.read(data_folder+'v4_'+dataset_long+'/glf_seed323_split1_nGPCgrid233.h5ad').var.index ),
              'genes236': np.array(sc.read(data_folder+'v4_'+dataset_long+'/glf_seed326_split1_nGPCgrid236.h5ad').var.index ),
              'genes239': np.array(sc.read(data_folder+'v4_'+dataset_long+'/glf_seed329_split1_nGPCgrid239.h5ad').var.index ) })

df_genes_ctrl.to_csv(data_folder+'v4_'+dataset_long+'/glf_nGPC_genes_ctrl.csv')

# len(np.intersect1d(genes227, genes230)) = 51
# len(np.intersect1d(genes227, genes233)) = 56
# len(np.intersect1d(genes227, genes236)) = 54
# len(np.intersect1d(genes227, genes239)) = 55


# len(np.intersect1d(genes230, genes233)) = 57
# len(np.intersect1d(genes230, genes236)) = 51
# len(np.intersect1d(genes230, genes239)) = 53

# len(np.intersect1d(genes233, genes236)) = 61
# len(np.intersect1d(genes233, genes239)) = 52

# len(np.intersect1d(genes236, genes239)) = 51
