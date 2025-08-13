dataset_long = 'greenleaf'
dataset_short = 'glf'
sct_seed = 615
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

gene_set_prefix = 'nGPCblk'

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *
from v4_functions import *

def compute_sct_avg_velocity(tnode,timesteps):
    v_shape = tnode.adata.shape
    v = np.zeros(v_shape)
    for t in timesteps:
        v += compute_sctour_velocity(tnode, timestep=t)
    return v/len(timesteps)

def sct_add_output(adata, tnode, timesteps):
    adata.layers['velocity'] = compute_sct_avg_velocity(tnode, timesteps)
    get_umap_sct(adata)

def sct_add_velo(split_seed, grid_seed, sct_seed=615):
    gene_set_name = gene_set_prefix+str(grid_seed)
    method = 'sct_'+gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
    print('############# read data')
    adata_prefix = 'adata_'+dataset_short+'_'+method
    tnode_prefix = 'tnode_'+dataset_short+'_'+method
    total = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_'+'total'+'_v4.h5ad') # 
    split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_'+'split1'+'_v4.h5ad') # 
    split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_'+'split2'+'_v4.h5ad') # 
    tnode_total = sct.predict.load_model(data_folder+'tnode_'+dataset_short+'_'+method+'_total_v4.pth')
    tnode_split1 = sct.predict.load_model(data_folder+'tnode_'+dataset_short+'_'+method+'_split1_v4.pth')
    tnode_split2 = sct.predict.load_model(data_folder+'tnode_'+dataset_short+'_'+method+'_split2_v4.pth')
    # average velocity
    print('############# compute average velocity')
    timesteps=[i/50 for i in range(1,11)]
    np.random.seed(sct_seed)
    random.seed(sct_seed)
    sct_add_output(total, tnode_total, timesteps)
    sct_add_output(split1, tnode_split1, timesteps)
    sct_add_output(split2, tnode_split2, timesteps)
    total.write(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
    split1.write(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
    split2.write(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

def sct_create_plots(split_seed, grid_seed, sct_seed=615):
    gene_set_name = gene_set_prefix+str(grid_seed)
    method = 'sct_'+gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
    adata_prefix = 'adata_'+dataset_short+'_'+method
    tnode_prefix = 'tnode_'+dataset_short+'_'+method
    print('############# read data')
    total = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_'+'total'+'_v4_outputAdded.h5ad') # 
    total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf']
    split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_'+'split1'+'_v4_outputAdded.h5ad') # 
    split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_'+'split2'+'_v4_outputAdded.h5ad') # 
    print('############# vector field')
    plot_vf_umap(adata_in=split1, data_version="split1",data=dataset_short,method=method,fig_folder=fig_folder)
    plot_vf_umap(adata_in=split2, data_version="split2",data=dataset_short,method=method,fig_folder=fig_folder)
    plot_vf_umap(adata_in=total, data_version="total",data=dataset_short,method=method,fig_folder=fig_folder)
    print('############# velocity')
    plot_sct_velocity(adata_in=split1,data_version='split1',dataset=dataset_short,fig_folder=fig_folder,method=method)
    plot_sct_velocity(adata_in=split2,data_version='split2',dataset=dataset_short,fig_folder=fig_folder,method=method)
    plot_sct_velocity(adata_in=total,data_version='total',dataset=dataset_short,fig_folder=fig_folder,method=method)
    print('############# cosine similarity')
    #c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
    #c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    print('############# confidence')
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)
    print('############# all done')

"""for i in range(5):
    split_seed = [317,320,323,326,329][i]
    grid_seed = [227,230,233,236,239][i]
    sct_create_plots(split_seed=split_seed, grid_seed=grid_seed)
    print('################## seed'+str(split_seed)+' done')
"""
"""
for i in range(4):
    split_seed = [320,323,326,329][i]
    grid_seed = [230,233,236,239][i]
    sct_add_velo(split_seed=split_seed, grid_seed=grid_seed)
    print('################## seed'+str(split_seed)+' done')
"""

for i in range(5):
    split_seed = [317,320,323,326,329][i]
    grid_seed = [227,230,233,236,239][i]
    sct_create_plots(split_seed=split_seed, grid_seed=grid_seed)
    print('################## seed'+str(split_seed)+' done')


print('################# all done')
