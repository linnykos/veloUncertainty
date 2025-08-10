import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

import numpy as np
import matplotlib.pyplot as plt

def compute_velocity_union_df(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    velo_split1 = pd.DataFrame(adata_split1.layers['velocity'], columns=velo_genes_split1)
    velo_split2 = pd.DataFrame(adata_split2.layers['velocity'], columns=velo_genes_split2)
    if method=='scv':
        velo_genes_split1 = velo_genes_split1[~np.isnan(velo_split1.loc[0])] 
        velo_genes_split2 = velo_genes_split2[~np.isnan(velo_split2.loc[0])] 
    union_genes_velo = np.union1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Size of the union of genes for velocity computation in splits = '+str(union_genes_velo.shape[0])) 
    Nrow = adata_split1.shape[0]
    velo_df1 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split1:
        velo_df1[gene] = velo_split1[gene]
    velo_df2 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split2:
        velo_df2[gene] = velo_split2[gene]
    return velo_df1, velo_df2

# velovi_woprep
split_seed = 317
method = 'velovi_woprep'
dataset_short = 'glf'
dataset_long = 'greenleaf'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)

df1, df2 = compute_velocity_union_df(split1,split2,method=method)

from sklearn.metrics.pairwise import cosine_similarity
cos_sim = cosine_similarity(df1,df2)

import numpy as np
np.where(split1.obsp['connectivities'][0].todense()!=0)[1]
idx_cell = 0

def compute_cos_sim_with_nb(cos_sim, split1):
    mean_nb = np.zeros(cos_sim.shape[1])
    mean_nnb = np.zeros(cos_sim.shape[1])
    for idx_cell in np.arange(cos_sim.shape[1]):
        idx_nb = np.where(split1.obsp['connectivities'][idx_cell].todense()!=0)[1]
        mean_nb[idx_cell] = np.mean( cos_sim[idx_cell,idx_nb] )
        mean_nnb[idx_cell] = np.mean( cos_sim[idx_cell, np.isin(np.arange(cos_sim.shape[1]), idx_nb, invert=True) ] )


