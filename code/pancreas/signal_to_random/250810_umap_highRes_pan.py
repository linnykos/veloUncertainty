# plot pancreas, high resolution
def plot_highRes_uamp_pamcreas(method, outputAdded):
    # setup
    dataset_long = 'pancreas'
    dataset_short = 'pan'
    split_seed=317
    celltype_label = get_celltype_label(dataset_short)
    data_version = 'total'
    # read data
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    data_method = dataset_short+"_"+method
    total_path = ''
    if outputAdded:
        total_path = data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad'
    else:
        total_path = data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4.h5ad'
    print(total_path)
    total = sc.read_h5ad(total_path)
    total.obsm['X_umap'] = total.obsm['X_umapOriginal'].copy()
    scv.pl.velocity_embedding_stream(total, basis='umap',color=celltype_label,recompute=True, legend_loc='none', arrow_size=2,
                                 title=dataset_short+'+'+method, dpi=300, 
                                 save=fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_highRes.png")    
    print(fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_highRes.png")

import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *



plot_highRes_uamp_pamcreas(method='scv', outputAdded=False)
plot_highRes_uamp_pamcreas(method='utv', outputAdded=False)
plot_highRes_uamp_pamcreas(method='sct', outputAdded=True)
plot_highRes_uamp_pamcreas(method='velovi', outputAdded=True)
plot_highRes_uamp_pamcreas(method='velovi_woprep', outputAdded=True)



def plot_highRes_uamp_pamcreasINC(method, outputAdded):
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    data_method = dataset_short+"_"+method
    total_path = ''
    if outputAdded:
        total_path = data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4_outputAdded.h5ad'
    else:
        total_path = data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4.h5ad'
    print(total_path)
    total = sc.read_h5ad(total_path)
    total.obsm['X_umap'] = total.obsm['X_umapOriginal'].copy()
    scv.pl.velocity_embedding_stream(total, basis='umap',color=celltype_label,recompute=True, legend_loc='none', arrow_size=2,
                                 title=dataset_short+'+'+method, dpi=300, 
                                 save=fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_highRes.png")    
    print(fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_highRes.png")



