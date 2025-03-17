method_prefix = 'velovi_woprep'
gene_set_name = 'Mark'
split_seed = 317

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


def run_velovi_ery_Mark(data_version, method_prefix, gene_set_name, split_seed):
    split_version = data_version
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method = method_prefix + '_' + gene_set_name
    #celltype_label = 'celltype'
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'
    from velovi import preprocess_data, VELOVI
    adata = None
    print_message_with_time("#################### "+dataset_long+'_'+method+'_'+data_version+": Read data ")
    if split_version=='total':
        adata = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_total_allgenes.h5ad')
    elif 'split' in split_version:
        adata = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_'+split_version+'_allgenes.h5ad')
    print_message_with_time("#################### "+data_version+": Preprocess data ")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if (not 'wop' in method): adata = preprocess_data(adata)
    # train and apply model
    print_message_with_time("#################### "+data_version+": Train model ")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    # save vae
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    print_message_with_time("#################### "+data_version+": Save vae ")
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'.pt',overwrite=True)
    # write data
    print_message_with_time("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'.h5ad')
    print_message_with_time("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)


def plot_velovi_ery_Mark(method_prefix, gene_set_name, split_seed, plot_total=True):
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method = method_prefix + '_' + gene_set_name
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    ## read data
    split1 = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split1'+'.h5ad')
    split2 = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split2'+'.h5ad')
    total = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'total'+'.h5ad')
    vae_split1 = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'split1.pt', split1)
    vae_split2 = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'split2.pt', split2)
    vae_total = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'total.pt', total)
    ## add velovi outputs to adata
    print_message_with_time("############## Add velovi outputs to adata")
    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    add_velovi_outputs_to_adata(total, vae_total)
    ######################################################
    ## compute umap
    print_message_with_time("############## Compute umap")
    compute_umap_ery(total)
    compute_umap_ery(split1)
    compute_umap_ery(split2)
    #split1.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split1'+'_outputAdded.h5ad')
    #split2.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split2'+'_outputAdded.h5ad')
    #total.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'total'+'_outputAdded.h5ad')
    ######################################################
    ## plot velocity
    scv.tl.velocity_graph(total)
    scv.tl.velocity_graph(split1)
    scv.tl.velocity_graph(split2)
    if (plot_total): plot_velocity(adata_in=total,fig_folder=fig_folder,data_version="total",dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split1,fig_folder=fig_folder,data_version="split1",dataset=dataset_short,method=method,split_seed=split_seed)
    plot_velocity(adata_in=split2,fig_folder=fig_folder,data_version="split2",dataset=dataset_short,method=method,split_seed=split_seed)
    print_message_with_time("############## plot velocity")
    ######################################################
    ## plot cosine similarity
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    print_message_with_time("############## plot cosine similarity")
    ## plot velo_conf
    scv.tl.velocity_confidence(total)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)
    print_message_with_time("############## plot velocity confidence")

"""
run_velovi_ery_Mark(data_version='total', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
run_velovi_ery_Mark(data_version='split1', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
run_velovi_ery_Mark(data_version='split2', method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)
"""
plot_velovi_ery_Mark(method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed)

