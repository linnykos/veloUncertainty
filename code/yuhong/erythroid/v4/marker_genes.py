"""
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

adata = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad')
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)
scv.pp.normalize_per_cell(adata) # normalize
scv.pp.log1p(adata) # log
sc.tl.leiden(adata)

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.settings.figdir = fig_folder
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='marker_gene.png')
adata.write(data_folder+'ery_total_leiden_wilcoxon.h5ad')

result = adata.uns['rank_genes_groups']

names = np.array(result['names'].tolist())
logfoldchanges = np.array(result['logfoldchanges'].tolist())
pvals_adj = np.array(result['pvals_adj'].tolist())
# Create a DataFrame for further filtering
markers = pd.DataFrame({
    'gene': names.flatten(),
    'logfoldchange': logfoldchanges.flatten(),
    'pvals_adj': pvals_adj.flatten()
})

markers.to_csv(data_folder+'ery_rank_genes_groups.csv')

dataset_long = 'erythroid'
dataset_short = 'ery'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

markers = pd.read_csv(data_folder+'ery_rank_genes_groups.csv')
marker_genes = markers[
    (markers['logfoldchange'].abs() >= 3) &  # |log2FC| ≥ 1
    (markers['pvals_adj'] <= 0.01)        # Adjusted p-value ≤ 0.05
]
"""
###################################################
## Mar.8
import numpy as np
import pandas as pd
import scanpy as sc

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

#adata = sc.read_h5ad(data_folder+dataset_short+'_total_allgenes.h5ad')
adata = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad')
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)
scv.pp.normalize_per_cell(adata) # normalize
scv.pp.log1p(adata) # log

sc.tl.rank_genes_groups(adata, 'celltype', method='wilcoxon')
sc.settings.figdir = fig_folder
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='marker_gene.png')
adata.write(data_folder+'ery_total_leiden_wilcoxon.h5ad')

result = adata.uns['rank_genes_groups']
names = np.array(result['names'].tolist())
logfoldchanges = np.array(result['logfoldchanges'].tolist())
pvals_adj = np.array(result['pvals_adj'].tolist())
# Create a DataFrame for further filtering
markers = pd.DataFrame({
    'gene': names.flatten(),
    'logfoldchange': logfoldchanges.flatten(),
    'pvals_adj': pvals_adj.flatten()
})

markers.to_csv(data_folder+'ery_rank_genes_groups.csv') # not needed


import collections
gene_counts = collections.Counter(markers['gene'][0:1000])
genes_selected = {gene for gene, count in gene_counts.items() if count > 1} # 265 genes
df_genes_selected = pd.DataFrame( { 'gene': np.array(list(genes_selected), dtype=object) } ) # convert into data frame

df_genes_selected.to_csv(data_folder+'ery_marker_genes_wilcoxon.csv')


###########################################################
## create data
import numpy as np
import pandas as pd
import scanpy as sc

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

gene_set_name = 'Mark'

total = sc.read_h5ad(data_folder+dataset_short+'_total_allgenes.h5ad') # total adata
df_genes_selected = pd.read_csv(data_folder+'ery_marker_genes_wilcoxon.csv')['gene'] # the selected marker genes
idx_marker_genes = [np.where(total.var.index==x)[0][0] for x in df_genes_selected] # index in adata objects
total[:,idx_marker_genes].write(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_total_allgenes.h5ad')

for split_seed in [317,320,323,326,329]:
    adata_split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
    adata_split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')
    adata_split1[:,idx_marker_genes].write(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_split1_allgenes.h5ad')
    adata_split2[:,idx_marker_genes].write(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_split2_allgenes.h5ad')
    # for utv
    adata_split1.var['highly_variable'] = True
    adata_split2.var['highly_variable'] = True
    adata_split1[:,idx_marker_genes].write(data_folder+'utv0_adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_split1_allgenes.h5ad')
    adata_split2[:,idx_marker_genes].write(data_folder+'utv0_adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_split2_allgenes.h5ad')
    print('seed'+str(split_seed)+' done!')

total = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_total_allgenes.h5ad')
total.var['highly_variable'] = True
total.write(data_folder+'utv0_adata_'+dataset_short+'_'+gene_set_name+'_total_allgenes.h5ad')

# check the subset of adata
# np.intersect1d(total[:,idx_marker_genes].var.index, df_genes_selected)
# np.intersect1d(adata_split1[:,idx_marker_genes].var.index, df_genes_selected).shape
# np.intersect1d(adata_split2[:,idx_marker_genes].var.index, df_genes_selected).shape

###########################################
###########################################
## nMark227
import scanpy as sc
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors
import pandas as pd

## functions
def select_genes_ctrl(df_tgt, df_ctrl, n_sample_scale, grid_seed, bs=4, w=1.0):
    # assign memberships to genes
    bins1 = np.linspace(np.min(df_tgt['fnzS']), np.max(df_ctrl['fnzS']), bs) # Step 1: Generate bins
    bins2 = np.linspace(np.min(df_tgt['fnzU']), np.max(df_ctrl['fnzU']), bs)
    bins3 = np.linspace(np.min(df_tgt['r2']), np.max(df_ctrl['r2']), bs)
    bin_data = np.array([[b1, b2, b3] for b1 in bins1 for b2 in bins2 for b3 in bins3]) 
    bin_dist = cdist(bin_data, bin_data, metric='euclidean') # Step 2: Compute Euclidean distances / Pairwise distances
    nbrs = NearestNeighbors(n_neighbors=1).fit(bin_data) # Step 4: Compute bin membership using kNN (k=1)
    membership_tgt = nbrs.kneighbors(df_tgt[['fnzS','fnzU','r2']], return_distance=False).flatten()
    membership_ctrl = nbrs.kneighbors(df_ctrl[['fnzS','fnzU','r2']], return_distance=False).flatten()
    unique, counts = np.unique(membership_tgt, return_counts=True) # the frequency table of target genes' membership
    freq_tgt = dict(zip(unique, counts))
    # clustering
    np.random.seed(grid_seed)
    genes_select = np.array([])
    for lab in freq_tgt.keys():
        print(lab)
        n_sample = freq_tgt[lab]*n_sample_scale
        idx_ctrl = (np.where( membership_ctrl==lab )[0])
        idx_ctrl = idx_ctrl[~np.isin(idx_ctrl, genes_select)]
        if (len(idx_ctrl) <= n_sample):
            genes_select = np.append( genes_select, idx_ctrl ) 
            n_sample -= len(idx_ctrl)
            nbr_list = np.array([lab]) 
            while (n_sample>0):
                print('n_sample='+str(n_sample))
                val = bin_dist[lab,]
                val[nbr_list] = np.Inf
                lab_nbr = (np.random.choice( np.where(val==np.min(val))[0], size=1 )[0])
                nbr_list = np.append(nbr_list, lab_nbr)
                idx_ctrl = (np.where( membership_ctrl==lab_nbr )[0])
                idx_ctrl = idx_ctrl[~np.isin(idx_ctrl, genes_select)]
                if len(idx_ctrl)==0:
                    print('->' + str(lab_nbr)+', case 1')
                    continue
                elif len(idx_ctrl) <= n_sample:
                    print('->' + str(lab_nbr)+', case 2')
                    genes_select = np.append(genes_select, idx_ctrl)
                    n_sample -= len(idx_ctrl)
                else:
                    print('->' + str(lab_nbr)+', case 3, len(idx_ctrl)>n_sample:', str(len(idx_ctrl)>n_sample))
                    idx_sample = np.random.choice(idx_ctrl, size=n_sample, replace=False, p=None)
                    genes_select = np.append(genes_select, idx_sample)
                    n_sample = 0
        else:
            idx_sample = np.random.choice(idx_ctrl, size=n_sample, replace=False, p=None)
            genes_select = np.append(genes_select, idx_sample)
    return genes_select

def get_adata_ctrl(df_tgt, df_ctrl, adata_ctrl, split_seed, grid_seed, n_sample_scale, bs=4, w=1.0):
    genes_ctrl_selected = select_genes_ctrl(df_tgt, df_ctrl, n_sample_scale, grid_seed, bs=bs, w=w)
    adata_ctrl_subset = adata_ctrl[:, [int(x) for x in genes_ctrl_selected]]
    split1 = sc.read(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
    split2 = sc.read(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')
    split1_ctrl = split1[:, adata_ctrl.var.index]
    split2_ctrl = split2[:, adata_ctrl.var.index]
    split1_ctrl_subset = split1_ctrl[:, [int(x) for x in genes_ctrl_selected]]
    split2_ctrl_subset = split2_ctrl[:, [int(x) for x in genes_ctrl_selected]]
    adata_ctrl_subset.var['highly_variable'] = True
    split1_ctrl_subset.var['highly_variable'] = True
    split2_ctrl_subset.var['highly_variable'] = True
    return adata_ctrl_subset, split1_ctrl_subset, split2_ctrl_subset


# 0. setup
dataset_long = 'erythroid'
dataset_short = 'ery'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
#fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'
#df_genes_selected = pd.read_csv(data_folder+'ery_marker_genes_wilcoxon.csv')
gene_mark = pd.read_csv(data_folder+'ery_marker_genes_wilcoxon.csv')['gene']

adata = sc.read_h5ad(data_folder+dataset_short+'_total_allgenes.h5ad') # total adata
genes = adata.var.index
total = adata.copy()

##################
## R2
indicator_mark = genes.isin( gene_mark )
## The likelihood computation requires recover_dynamics and is not available in stochastic or deterministic mode.
## to compute gene R2, it is unavoidable to filter out a few genes
fnzS_all = np.asarray(np.mean(adata.layers["spliced"] > 0, axis=0))[0]
fnzU_all = np.asarray(np.mean(adata.layers["unspliced"] > 0, axis=0))[0]

import scvelo as scv
## r2
scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode='deterministic')
r2_all = np.array(adata.var.velocity_r2)
df = pd.DataFrame( {'gene_name': genes, 'mark': indicator_mark, 'fnzS': fnzS_all, 'fnzU': fnzU_all, 'r2': r2_all} )

df.to_csv(data_folder+dataset_short+'_r2_df_mark.csv')
## r2 file saved
##################

df = pd.read_csv(data_folder+dataset_short+'_r2_df_mark.csv')
df_mark = df[df['mark']]
df_nmark = df[~df['mark']]

n_sample_scale = 2 # used to be 1
gene_set_name = 'nMark'
for i in range(5):
    split_seed = [317,320,323,326,329][i]
    grid_seed = [227,230,233,236,239][i]
    total_nMark, split1_nMark, split2_nMark = get_adata_ctrl(df_tgt=df_mark, df_ctrl=df_nmark,
                                                             adata_ctrl=total[:, df[~df['mark']]['gene_name']], 
                                                             split_seed=split_seed, grid_seed=grid_seed, n_sample_scale=n_sample_scale)
    total_nMark.write(data_folder+'adata_'+dataset_short+'_'+gene_set_name+str(grid_seed)+'_total.h5ad')
    split1_nMark.write(data_folder+'adata_'+dataset_short+'_'+gene_set_name+str(grid_seed)+'_split1.h5ad')
    split2_nMark.write(data_folder+'adata_'+dataset_short+'_'+gene_set_name+str(grid_seed)+'_split2.h5ad')
    print('################ '+str(split_seed)+' done')

# save the control genes
df_genes_ctrl = pd.DataFrame({'genes227': np.array(sc.read(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'227_split1.h5ad').var.index),
              'genes230': np.array(sc.read(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'230_split1.h5ad').var.index ),
              'genes233': np.array(sc.read(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'233_split1.h5ad').var.index ),
              'genes236': np.array(sc.read(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'236_split1.h5ad').var.index ),
              'genes239': np.array(sc.read(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'239_split1.h5ad').var.index ) })

df_genes_ctrl.to_csv(data_folder+dataset_short+'_'+gene_set_name+'_genes_ctrl.csv')

"""
np.intersect1d(df_genes_ctrl['genes227'],df_genes_ctrl['genes230']).shape # 222
np.intersect1d(df_genes_ctrl['genes227'],df_genes_ctrl['genes233']).shape # 217
np.intersect1d(df_genes_ctrl['genes227'],df_genes_ctrl['genes236']).shape # 228
np.intersect1d(df_genes_ctrl['genes227'],df_genes_ctrl['genes239']).shape # 213

np.intersect1d(df_genes_ctrl['genes230'],df_genes_ctrl['genes233']).shape # 220
np.intersect1d(df_genes_ctrl['genes230'],df_genes_ctrl['genes236']).shape # 226
np.intersect1d(df_genes_ctrl['genes230'],df_genes_ctrl['genes239']).shape # 223

np.intersect1d(df_genes_ctrl['genes233'],df_genes_ctrl['genes236']).shape # 221
np.intersect1d(df_genes_ctrl['genes233'],df_genes_ctrl['genes239']).shape # 226

np.intersect1d(df_genes_ctrl['genes236'],df_genes_ctrl['genes239']).shape # 226
"""
