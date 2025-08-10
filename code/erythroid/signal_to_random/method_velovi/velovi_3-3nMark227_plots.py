method_prefix = 'velovi'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_velovi import *
from v4_functions import *

def plot_velovi_ery_Mark(method_prefix, gene_set_name, split_seed, plot_total=True):
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method = method_prefix + '_' + gene_set_name
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    ## read data
    split1 = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split1'+'_outputAdded.h5ad')
    split2 = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split2'+'_outputAdded.h5ad')
    total = sc.read_h5ad(data_folder+'seed317'+'/'+method+'/'+'adata_'+dataset_short+'_'+method+'_'+'total'+'_outputAdded.h5ad')
    #vae_split1 = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'split1.pt', split1)
    #vae_split2 = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'split2.pt', split2)
    #vae_total = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'total.pt', total)
    ## add velovi outputs to adata
    ######################################################
    ## plot velocity
    print_message_with_time("############## plot velocity")
    if (plot_total): plot_velocity(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
    #plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
    #plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)
    ######################################################
    ## plot cosine similarity
    print_message_with_time("############## plot cosine similarity")
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,
                           fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,
                                            method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,
                                               method=method,fig_folder=fig_folder,split_seed=split_seed)
    ## plot velo_conf
    print_message_with_time("############## plot velocity confidence")
    scv.tl.velocity_confidence(total)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)
    print_message_with_time("############## all done")

grid_seed = 227
gene_set_name = 'nMark' + str(grid_seed)
for i in range(4):
    split_seed = [320, 323, 326, 329][i]
    plot_velovi_ery_Mark(method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed, plot_total=False)

