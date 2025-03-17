import sctour as sct
import scanpy as sc
import numpy as np
import torch
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *
from v4_functions import print_message_with_time

def run_sct_ery_Mark(gene_set_name,split_version,split_seed,sct_seed=615):
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method_prefix = 'sct'
    method = method_prefix + '_' + gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    print_message_with_time("################## Read data")
    raw = read_raw_adata(dataset_short)
    if split_version=='total':
        adata = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_total_allgenes.h5ad')
    elif 'split' in split_version:
        adata = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+gene_set_name+'_seed'+str(split_seed)+'_'+split_version+'_allgenes.h5ad')
    adata.obsm['X_umapOriginal'] = raw.obsm['X_umap'].copy()
    adata.obsm['X_pcaOriginal'] = raw.obsm['X_pca'].copy()
    print_message_with_time("########### start to train model for "+split_version+' ')
    tnode = sct_train_and_return_tnode(adata, sct_seed)
    print_message_with_time("########### start to compute velocity for "+split_version+' ')
    print_message_with_time("########### "+split_version+' velocity computed, start to write data')
    adata.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+split_version+'.h5ad')
    tnode.save_model(save_dir=savedata_folder, save_prefix='tnode_'+dataset_short+'_'+method+'_'+split_version)
    print_message_with_time("########### "+split_version+' data wrote')

##########################################
## plot
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from sctour_misc import *


def compute_sct_avg_velocity(tnode,timesteps):
    v_shape = tnode.adata.shape
    v = np.zeros(v_shape)
    for t in timesteps:
        v += compute_sctour_velocity(tnode, timestep=t)
    return v/len(timesteps)

def plot_sct_ery_Mark(gene_set_name, split_seed, plot_total=True):
    import scvelo as scv
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method_prefix = 'sct'
    method = method_prefix + '_' + gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
    tnode_total = sct.predict.load_model(data_folder+'seed'+str(split_seed)+'/'+method+'/tnode_'+dataset_short+'_'+method+'_total.pth')
    tnode_split1 = sct.predict.load_model(data_folder+'seed'+str(split_seed)+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split1.pth')
    tnode_split2 = sct.predict.load_model(data_folder+'seed'+str(split_seed)+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split2.pth')
    # recover velocity
    timesteps=[i/50 for i in range(1,11)]
    total.layers['velocity'] = compute_sct_avg_velocity(tnode_total, timesteps)
    split1.layers['velocity'] = compute_sct_avg_velocity(tnode_split1, timesteps) 
    split2.layers['velocity'] = compute_sct_avg_velocity(tnode_split2, timesteps)
    print('Velocity computed')
    # compute umap
    get_umap_sct(total)
    get_umap_sct(split1)
    get_umap_sct(split2)
    print('UMAP computed')
    # vector field
    plot_vf_umap(adata_in=split1, data_version="split1",data=dataset_short,method=method,fig_folder=fig_folder)
    plot_vf_umap(adata_in=split2, data_version="split2",data=dataset_short,method=method,fig_folder=fig_folder)
    if (plot_total): plot_vf_umap(adata_in=total, data_version="total",data=dataset_short,method=method,fig_folder=fig_folder)
    # velocity
    plot_sct_velocity(adata_in=split1,data_version='split1',dataset=dataset_short,fig_folder=fig_folder)
    plot_sct_velocity(adata_in=split2,data_version='split2',dataset=dataset_short,fig_folder=fig_folder)
    if (plot_total):
        plot_sct_velocity(adata_in=total,data_version='total',dataset=dataset_short,fig_folder=fig_folder)
    ## plot cosine similarity
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,
                           fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,
                                            method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,
                                               method=method,fig_folder=fig_folder,split_seed=split_seed)
    ######################################################
    ## plot velo_conf
    scv.tl.velocity_confidence(total)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)

"""
run_sct_ery_Mark(gene_set_name='Mark', split_version='total', split_seed=317,sct_seed=615)
run_sct_ery_Mark(gene_set_name='Mark', split_version='split1', split_seed=317,sct_seed=615)
run_sct_ery_Mark(gene_set_name='Mark', split_version='split2', split_seed=317,sct_seed=615)
"""
plot_sct_ery_Mark(gene_set_name='Mark', split_seed=317, plot_total=True)


