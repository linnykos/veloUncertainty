import scvelo as scv
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *


dataset_long='pancreas'
dataset_short='pan'
method='velovi_woprep'
split_seed=317
total = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='total',allgenes=False,outputAdded=True)
np.mean(total.layers['velocity'], axis=0)
# axis=0 mean over row, axis=1 mean over column

adata = total

celltype_label = 'clusters'

# deviation to the mean within one velocity set
for ct in adata.obs[celltype_label].cat.categories:
    ct_idx = np.where(adata.obs[celltype_label]==ct)[0]
    ct_mean_velo = np.mean(adata.layers['velocity'][ct_idx,:], axis=0)
    ct_cos_sim_with_mean = np.diag(cosine_similarity(adata.layers['velocity'][ct_idx,:], np.tile(ct_mean_velo, (len(ct_idx), 1))))
    print(ct)
    print( 'mean:' + str(np.round(np.mean(ct_cos_sim_with_mean),4)) + ', median:'+ str(np.round(np.median(ct_cos_sim_with_mean),4)))

dataset_long='pancreas'
dataset_short='pan'
celltype_label = 'clusters'
method='utv'
split_seed=317


s1 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split1',allgenes=False,outputAdded=True)
s2 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split2',allgenes=False,outputAdded=True)

def compute_cos_sim_rm_celltype_mean(s1,s2):
    # between splits
    celltypes = s1.obs[celltype_label].cat.categories
    for ct in celltypes:
        ct_idx = np.where(s1.obs[celltype_label]==ct)[0]
        print(ct+', N='+str(len(ct_idx)))
        ct_velo_s1 = s1.layers['velocity'][ct_idx,:]
        ct_velo_s2 = s2.layers['velocity'][ct_idx,:]
        ct_velo_mean_s1 = np.mean(ct_velo_s1, axis=0)
        ct_velo_mean_s2 = np.mean(ct_velo_s2, axis=0)
        # cannot do this since genes are not matched
        #ct_mean_velo_cos_sim = np.diag(cosine_similarity(ct_mean_velo_s1.reshape(1,-1),ct_mean_velo_s2.reshape(1,-1)))[0]
        ct_velo_rm_mean_s1 = ct_velo_s1 - ct_velo_mean_s1
        ct_velo_rm_mean_s2 = ct_velo_s2 - ct_velo_mean_s2
        velo_genes_s1 = s1.var.index
        velo_genes_s2 = s2.var.index
        velo_s1 = pd.DataFrame(ct_velo_rm_mean_s1, columns=velo_genes_s1)
        velo_s2 = pd.DataFrame(ct_velo_rm_mean_s2, columns=velo_genes_s2)
        if method=='scv':
            velo_genes_s1 = velo_genes_s1[~np.isnan(velo_s1.loc[0])] 
            velo_genes_s2 = velo_genes_s2[~np.isnan(velo_s2.loc[0])] 
        union_genes_velo = np.union1d(np.array(velo_genes_s1), np.array(velo_genes_s2))
        Nrow = len(ct_idx)
        velo_df1 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
        for gene in velo_genes_s1: velo_df1[gene] = velo_s1[gene]
        velo_df2 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
        for gene in velo_genes_s2: velo_df2[gene] = velo_s2[gene]
        ct_cos_sim_rm_mean = np.diag(cosine_similarity(velo_df1,velo_df2))
        print( 'mean:' + str(np.round(np.mean(ct_cos_sim_rm_mean),4)) + ', median:'+ str(np.round(np.median(ct_cos_sim_rm_mean),4)))

print_velo_rm_mean_cos_sim_btw_splits(s1,s2)

s1_tmp = s1.copy()[np.where(s1.obs[celltype_label]=='Epsilon')[0],]
s2_tmp = s2.copy()[np.where(s1.obs[celltype_label]=='Epsilon')[0],]
s1_tmp.layers['velocity'] = s1_tmp.layers['velocity']-np.mean(s1_tmp.layers['velocity'], axis=0)
s2_tmp.layers['velocity'] = s2_tmp.layers['velocity']-np.mean(s2_tmp.layers['velocity'], axis=0)
np.quantile(compute_cosine_similarity_union(s1_tmp,s2_tmp,method)[0], [0.,.25,.5,.75,1.])

cosine_similarity(np.tile(np.mean(total.layers['velocity'], axis=0),(1,1)), np.tile(ct_mean_velo,(1,1)))


dataset_long='pancreas'
dataset_short='pan'
method='velovi_woprep'
split_seed=317

s1 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split1',allgenes=False,outputAdded=True)
s2 = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split2',allgenes=False,outputAdded=True)

for ct in s1.obs[celltype_label].cat.categories:
    print(ct+', N='+str(len(ct_idx)))
    ct_idx = np.where(s1.obs[celltype_label]==ct)[0]
    s1_tmp = s1.copy()[ct_idx,:]
    s2_tmp = s2.copy()[ct_idx,:]
    ct_velo_s1 = s1_tmp.layers['velocity']
    ct_velo_s2 = s2_tmp.layers['velocity']
    ct_velo_mean_s1 = np.mean(ct_velo_s1, axis=0)
    ct_velo_mean_s2 = np.mean(ct_velo_s2, axis=0)
    s1_tmp.layers['velocity'] = s1_tmp.layers['velocity']-ct_velo_mean_s1
    s2_tmp.layers['velocity'] = s2_tmp.layers['velocity']-ct_velo_mean_s2
    print(np.quantile(compute_cosine_similarity_union(s1_tmp,s2_tmp,method)[0], [0.,.25,.5,.75,1.]))

#########################
"""
# copy velocity from split_velo to split_neigh
def get_velocity_from_twin_split(split_velo,split_neigh,method):
    velo_genes_v = split_velo.var.index
    velo_genes_n = split_neigh.var.index
    velo_split_v = pd.DataFrame(split_velo.layers['velocity'], columns=velo_genes_v)
    #velo_split_n = pd.DataFrame(split_neigh.layers['velocity'], columns=velo_genes_n)
    if method=='scv':
        res = pd.DataFrame(np.full(split_neigh.shape, np.nan), columns=velo_genes_n)
    else: 
        res = pd.DataFrame(np.zeros(split_neigh.shape), columns=velo_genes_n)
    for gene in velo_genes_v: 
        if gene in res.columns:
            res[gene] = velo_split_v[gene]
    return np.array(res)
"""

def copy_neighbor_from_twin_split(split_to,split_from):
    split = split_to.copy()
    split.uns['pca'] = split_from.uns['pca'].copy()
    split.uns['umap'] = split_from.uns['umap']
    split.obsm['X_pca'] = split_from.obsm['X_pca']
    split.obsm['X_umap'] = split_from.obsm['X_umap']
    split.obsp['connectivities'] = split_from.obsp['connectivities'].copy()
    split.obsp['distances'] = split_from.obsp['distances'].copy()
    scv.tl.velocity_graph(split, n_jobs=8)
    return split

def compute_confidence_cross_split_neigh(dataset_long,dataset_short,method,split_seed,outputAdded=False):
    split1_safe = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split1',allgenes=False,outputAdded=outputAdded)
    split2_safe = read_data_v4(dataset_long=dataset_long,dataset_short=dataset_short,method=method,split_seed=split_seed,data_version='split2',allgenes=False,outputAdded=outputAdded)
    scv.tl.velocity_confidence(split1_safe)
    scv.tl.velocity_confidence(split2_safe)
    conf11 = split1_safe.obs['velocity_confidence'] # confidence: neigh from split1, velo from split1
    conf22 = split2_safe.obs['velocity_confidence'] # confidence: neigh from split2, velo from split2
    # copy split1 velocities into split2
    split1 = copy_neighbor_from_twin_split(split_to=split1_safe,split_from=split2_safe)
    split2 = copy_neighbor_from_twin_split(split_to=split2_safe,split_from=split1_safe)
    #scv.tl.velocity_graph(split1, n_jobs=8)
    scv.tl.velocity_confidence(split1)
    #scv.tl.velocity_graph(split2, n_jobs=8)
    scv.tl.velocity_confidence(split2)
    conf12 = split1.obs['velocity_confidence'] # confidence: neigh from split1, velo from split2 
    conf21 = split2.obs['velocity_confidence'] # confidence: neigh from split2, velo from split1
    return [conf11,conf12,conf21,conf22]

def print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method,split_seed,outputAdded=True):
    conf11,conf12,conf21,conf22 = compute_confidence_cross_split_neigh(dataset_long,dataset_short,method,split_seed,outputAdded=outputAdded)
    print(dataset_short+'+'+method+', seed='+str(split_seed))
    print('cor(conf12,conf21) and [mean(conf12) mean(conf21)]')
    print( np.round( np.corrcoef(conf12,conf21)[0,1], 4 ) )
    print( np.round( (np.mean(conf12),np.mean(conf21)), 4) )
    #print( np.round(np.corrcoef([conf11,conf12,conf21,conf22]),5) )


## erythroid
dataset_long='erythroid'
dataset_short='ery'
split_seed=317
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed,outputAdded=False)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='utv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='sct',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi_woprep',split_seed=split_seed)




## larryMult
dataset_long='larryMult'
dataset_short='larryMult'
split_seed=317
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed,outputAdded=False)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='utv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='sct',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi_woprep',split_seed=split_seed)

split_seed=320
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed,outputAdded=False)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='utv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='sct',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi_woprep',split_seed=split_seed)

split_seed=323
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed,outputAdded=False)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='utv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='sct',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi_woprep',split_seed=split_seed)

split_seed=326
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed,outputAdded=False)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='utv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='sct',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi_woprep',split_seed=split_seed)

split_seed=329
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed,outputAdded=False)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='utv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='sct',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi_woprep',split_seed=split_seed)


## pancreasINC
dataset_long='pancreasINC'
dataset_short='panINC'
split_seed=317
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='utv',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='sct',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi',split_seed=split_seed)
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='velovi_woprep',split_seed=split_seed)




dataset_long='pancreas'
dataset_short='pan'
split_seed=317
print_confidence_summary_cross_split_neigh(dataset_long,dataset_short,method='scv',split_seed=split_seed)






np.round( np.corrcoef(conf12,conf21)[0,1], 5 ) 
np.round( (np.mean(conf12),np.mean(conf21)), 5)
np.corrcoef([conf11,conf12,conf21,conf22])
np.quantile(conf11-conf22,[0.,.25,.5,.75,1.]) #[-0.3126018 , -0.02132162, -0.00293523,  0.00279859,  0.13608122]