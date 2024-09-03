method = 'scv'
dataset_long = 'erythroid'

dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'scv'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

import scvelo as scv
import scanpy as sc
import bbknn
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

celltype_label = get_celltype_label(dataset_short)

split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_split1_v4.h5ad')

split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')
total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_total_v4_outputAdded.h5ad')

### feed neighborhood information in split1 to split2_allgenes and compute velocities
### compute cosine similarity between total and split2 with split1 neighbors
split2.obsp['connectivities'] = split1.obsp['connectivities'].copy()
split2.obsp['distances'] = split1.obsp['distances'].copy()
split2.uns['neighbors'] = split1.uns['neighbors'].copy()

import bbknn
scv.pp.filter_and_normalize(split2, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(split2, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(split2, batch_key='sequencing.batch')
split2.X = split2.X.toarray()
bbknn.ridge_regression(split2, batch_key='sample', confounder_key='celltype')
sc.tl.pca(split2)
bbknn.bbknn(split2, batch_key='sequencing.batch')
print("Batch correction done!")
sc.pp.neighbors(split2, n_neighbors=30, n_pcs=40) 
sc.tl.umap(split2)
scv.tl.recover_dynamics(split2,n_jobs=8)
scv.tl.velocity(split2, mode="dynamical")
scv.tl.velocity_graph(split2,n_jobs=8)

scv.pl.velocity_embedding_stream(split2, basis='umap',color='celltype',recompute=True,
                                     title='Velocity '+dataset_short+'+'+method+' (split_seed='+str(split_seed)+')',
                                     save=fig_folder+"velocity/test_"+dataset_short+'_'+method+"_umapCompute.png")
split2.obsm['X_umap'] = split1.obsm['X_umapOriginal'].copy()
scv.pl.velocity_embedding_stream(split2, basis='umap',color='celltype',recompute=True,
                                    title='Velocity '+dataset_short+'+'+method+' (split_seed='+str(split_seed)+')',
                                    save=fig_folder+"velocity/test_"+dataset_short+'_'+method+"_umapOriginal.png")    

split2.layers['velocity']
total.layers['velocity']

ci,ni = compute_cosine_similarity_intersect(adata_split1=split2,adata_split2=total,method=method)

cu,nu = compute_cosine_similarity_union(adata_split1=split2,adata_split2=total,method=method)

np.quantile(ci,[0.,.25,.5,.75,1.])
np.quantile(cu,[0.,.25,.5,.75,1.])

"""
>>> ci,ni = compute_cosine_similarity_intersect(adata_split1=split2,adata_split2=total,method=method)
Number of overlapped genes for velocity computation in splits = 268
>>> cu,nu = compute_cosine_similarity_union(adata_split1=split2,adata_split2=total,method=method)
Size of the union of genes for velocity computation in splits = 643
>>> np.quantile(ci,[0.,.25,.5,.75,1.])
[0.,.25,.5,.75,1.])
array([-0.40152761,  0.49668513,  0.58227361,  0.64651495,  0.82505851])
>>> np.quantile(cu,[0.,.25,.5,.75,1.])
array([-0.30606965,  0.41350181,  0.49077196,  0.54986151,  0.71526196])
"""
