method = 'sct_GPC'
dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed317/'+method+'/'

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *
from v4_functions import *

adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method

total = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed317/'+method+'/'+adata_prefix+'_total_v4_outputAdded.h5ad') # 
tnode_total = sct.predict.load_model('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed317/'+method+'/'+tnode_prefix+'_total_v4.pth')
total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()

def compute_sct_avg_velocity(tnode,timesteps):
    v_shape = tnode.adata.shape
    v = np.zeros(v_shape)
    for t in timesteps:
        v += compute_sctour_velocity(tnode, timestep=t)
    return v/len(timesteps)

for split_seed in [320,323,326,329]:
    print('######################### '+str(split_seed)+' starts!')
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
    split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4.h5ad') # 
    split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4.h5ad') # 
    tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
    tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')
    timesteps=[i/50 for i in range(1,11)]
    split1.layers['velocity'] = compute_sct_avg_velocity(tnode_split1, timesteps) 
    split2.layers['velocity'] = compute_sct_avg_velocity(tnode_split2, timesteps)
    get_umap_sct(split1)
    get_umap_sct(split2)
    print('UMAP computed')
    split1.write_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
    split2.write_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
    print('wrote data_outputAdded')
    """
    split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
    split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
    """
    ############################################
    # vector field
    print('############# vector field')
    #plot_vf_umap(adata_in=split1, data_version="split1",data=dataset_short,method=method,fig_folder=fig_folder)
    #plot_vf_umap(adata_in=split2, data_version="split2",data=dataset_short,method=method,fig_folder=fig_folder)
    plot_vf_umap(adata_in=total, data_version="total",data=dataset_short,method=method,fig_folder=fig_folder)
    ptime_sct_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)
    ############################################
    # velocity
    print('############# velocity')
    plot_sct_velocity(adata_in=split1,data_version='split1',dataset=dataset_short,fig_folder=fig_folder,method=method)
    plot_sct_velocity(adata_in=split2,data_version='split2',dataset=dataset_short,fig_folder=fig_folder,method=method)
    ######################################################
    ## plot cosine similarity
    print('############# cosine similarity')
    #c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
    #c2,n2 = compute_cosine_similarity_union(split1,split2,method)
    plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
    plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    ######################################################
    ## plot velo_conf
    print('############# confidence')
    scv.tl.velocity_confidence(total)
    plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
    plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)
    print('######################### '+str(split_seed)+' done!')



