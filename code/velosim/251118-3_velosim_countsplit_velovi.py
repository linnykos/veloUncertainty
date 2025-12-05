######
# velosim
import pandas as pd
import numpy as np
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/velosim/251118_test/'
S = pd.read_csv(data_folder+'velosim_counts_s.csv')
U = pd.read_csv(data_folder+'velosim_counts_u.csv')
overdisp_S = np.array(pd.read_csv(data_folder+'velosim_overdisp_S.csv')['x'])
overdisp_U = np.array(pd.read_csv(data_folder+'velosim_overdisp_U.csv')['x'])

split_seed = 317

np.random.seed(split_seed)
s1, s2  = countsplit(S,overdisps=overdisp_S)
u1, u2  = countsplit(U,overdisps=overdisp_U)

# create adata objects
import anndata as ad
import scipy.sparse as sp

velo_true = pd.read_csv(data_folder+'velosim_velocity.csv')
# total
adata = ad.AnnData(X = sp.csr_matrix(S))
adata.layers["spliced"] = sp.csr_matrix(S)
adata.layers["unspliced"] = sp.csr_matrix(U)
adata.layers['velocity_true'] = velo_true
# split1
adata1 = ad.AnnData(X = sp.csr_matrix(s1))
adata1.layers["spliced"] = sp.csr_matrix(s1)
adata1.layers["unspliced"] = sp.csr_matrix(u1)
adata.layers['velocity_true'] = velo_true
# split2
adata2 = ad.AnnData(X = sp.csr_matrix(s2))
adata2.layers["spliced"] = sp.csr_matrix(s2)
adata2.layers["unspliced"] = sp.csr_matrix(u2)
adata.layers['velocity_true'] = velo_true


method = 'velovi'
import scvelo as scv
def run_velovi_velosim(adata, method):
    # velovi
    from velovi import preprocess_data, VELOVI
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) 
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if (not 'wop' in method): adata = preprocess_data(adata)
    # train and apply model
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    return adata,vae

adata, vae = run_velovi_velosim(adata, method)
adata1, vae1 = run_velovi_velosim(adata1, method)
adata2, vae2 = run_velovi_velosim(adata2, method)

from v4_functions_velovi import *
adata = vae.adata.copy()
adata1 = vae1.adata.copy()
adata2 = vae2.adata.copy()

add_velovi_outputs_to_adata(adata1, vae1)
add_velovi_outputs_to_adata(adata2, vae2)
add_velovi_outputs_to_adata(adata, vae)


import scanpy as sc
method = 'velovi'
dataset_short = 'velosim251118'
dataset_long = 'velosim251118'
#sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)
compute_umap(adata1, dataset_short)
compute_umap(adata2, dataset_short)
compute_umap(adata, dataset_short)
adata.obsm['X_umapOriginal'] = adata.obsm['X_umap'].copy()
adata1.obsm['X_umapOriginal'] = adata1.obsm['X_umap'].copy()
adata2.obsm['X_umapOriginal'] = adata2.obsm['X_umap'].copy()


# make plots
import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

method = 'velovi'
dataset_short = 'velosim251118'
dataset_long = 'velosim251118'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/velosim/251118/'+method+'/'
plot_velocity(adata_in=adata,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity(adata_in=adata1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
plot_velocity(adata_in=adata2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)

plot_cosine_similarity(adata_split1=adata1,adata_split2=adata2,adata_total=adata,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(adata1,adata2,adata,dataset_short,method,fig_folder,split_seed)

scv.tl.velocity_confidence(adata)
scv.tl.velocity_confidence(adata1)
scv.tl.velocity_confidence(adata1)

plot_veloConf_and_cosSim(adata,adata1,adata2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(adata,dataset_short,method,fig_folder,split_seed)

# compare computed velocity with true velocity
np.diag(cosine_similarity(adata.layers['velocity'],adata.layers['velocity_true']))

np.mean( np.diag(cosine_similarity(adata.layers['velocity'],adata.layers['velocity_true'])) ) # 0.0138
np.median( np.diag(cosine_similarity(adata.layers['velocity'],adata.layers['velocity_true'])) ) # 0.0092


