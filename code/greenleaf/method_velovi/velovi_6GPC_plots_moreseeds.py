import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_velovi import *


method = 'velovi_GPC'
dataset_short = 'glf'
dataset_long = 'greenleaf'

data_total_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed317/'+method+'/'

total = sc.read_h5ad(data_total_folder+'adata_glf_velovi_GPC_total_GPC_outputAdded.h5ad')
vae_total = VELOVI.load(data_total_folder+'vae_glf_velovi_GPC_total_GPC.pt', total)
total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()

def compute_umap_GPC(adata):
    scv.pp.moments(adata, n_pcs=5, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=5) # used to be 10
    sc.tl.umap(adata)
    #scv.tl.velocity_graph(adata)

for split_seed in [320, 323, 326, 329]:
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

    split1 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split1_GPC.h5ad') # 
    split2 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split2_GPC.h5ad') # 

    vae_split1 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split1_GPC.pt', split1)
    vae_split2 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split2_GPC.pt', split2)
    #######################################
    ## add velovi outputs to adata
    print_message_with_time("############## Add velovi outputs to adata")
    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    ######################################################
    ## compute umap
    print_message_with_time("############## Compute umap")
    compute_umap_GPC(split1)
    compute_umap_GPC(split2)
    ## write data
    split1.write_h5ad(data_folder+'adata_glf_velovi_GPC_split1_GPC_outputAdded.h5ad')
    split2.write_h5ad(data_folder+'adata_glf_velovi_GPC_split2_GPC_outputAdded.h5ad')
    """
    split1 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split1_GPC_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split2_GPC_outputAdded.h5ad')

    vae_split1 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split1_GPC.pt', split1)
    vae_split2 = VELOVI.load(data_folder+'vae_glf_velovi_GPC_split2_GPC.pt', split2)
    """
    ######################################################
    ## plot cosine similarity
    print('######## cosine similarity')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) 
    c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    print(np.quantile(c1,[0.,.25,.5,.75,1.]) )
    print(np.quantile(c2,[0.,.25,.5,.75,1.]) )






