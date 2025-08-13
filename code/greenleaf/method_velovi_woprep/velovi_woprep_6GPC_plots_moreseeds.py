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

method = 'velovi_woprep_GPC'
dataset_short = 'glf'
dataset_long = 'greenleaf'

total = sc.read_h5ad("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed317/'+method+'/adata_glf_velovi_woprep_GPC_total_GPC_outputAdded.h5ad') # 
vae_total = VELOVI.load("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed317/'+method+'/vae_glf_velovi_woprep_GPC_total_GPC.pt', total)
#total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()


for split_seed in [320,323,326,329]:
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

    split1 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split1_GPC.h5ad') # 
    split2 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split2_GPC.h5ad') # 
    vae_split1 = VELOVI.load(data_folder+'vae_glf_'+method+'_split1_GPC.pt', split1)
    vae_split2 = VELOVI.load(data_folder+'vae_glf_'+method+'_split2_GPC.pt', split2)
    #######################################
    ## add velovi outputs to adata
    print_message_with_time("############## Add velovi outputs to adata")

    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    ######################################################
    ## compute umap
    print_message_with_time("############## Compute umap")
    compute_umap(split1, dataset_short)
    compute_umap(split2, dataset_short)
    ## write data
    split1.write_h5ad(data_folder+'adata_glf_'+method+'_split1_GPC_outputAdded.h5ad')
    split2.write_h5ad(data_folder+'adata_glf_'+method+'_split2_GPC_outputAdded.h5ad')
    """
    split1 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split1_GPC_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'adata_glf_'+method+'_split2_GPC_outputAdded.h5ad')
    vae_split1 = VELOVI.load(data_folder+'vae_glf_'+method+'_split1_GPC.pt', split1)
    vae_split2 = VELOVI.load(data_folder+'vae_glf_'+method+'_split2_GPC.pt', split2)
    """
    ######################################################
    ## plot velocity
    print('######## velocity')
    plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)

    ######################################################
    ## plot cosine similarity
    print('######## cosine similarity')
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)





