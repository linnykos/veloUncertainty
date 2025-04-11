dataset_long = 'pancreas'
dataset_short = 'pan'
gene_set_prefix = 'Mark'

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

def sct_create_plots(split_seed, grid_seed, sct_seed=615):
    gene_set_name = gene_set_prefix+str(grid_seed)
    method = 'sct_'+gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
    print('############# read data')
    adata_prefix = 'adata_'+dataset_short+'_'+method
    tnode_prefix = 'tnode_'+dataset_short+'_'+method
    total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4.h5ad') # 
    split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4.h5ad') # 
    split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4.h5ad') # 
    tnode_total = sct.predict.load_model(data_folder+tnode_prefix+'_total_v4.pth')
    tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
    tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')
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




"""
total.write_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1.write_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2.write_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
print('wrote data_outputAdded')


total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
"""

sct_create_plots(split_seed=317, grid_seed=227)
print('################## seed317 done')
sct_create_plots(split_seed=320, grid_seed=230)
print('################## seed320 done')
sct_create_plots(split_seed=323, grid_seed=233)
print('################## seed323 done')
sct_create_plots(split_seed=326, grid_seed=236)
print('################## seed326 done')
sct_create_plots(split_seed=329, grid_seed=239)
print('################## seed329 done')

"""
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
def sct_create_cos_sim_df_across_seed(velo_method, split_seeds, grid_seeds):
    df = pd.DataFrame()
    for i in range(len(split_seeds)):
        split_seed = split_seeds[i]
        grid_seed = grid_seeds[i]
        gene_set_name = gene_set_prefix+str(grid_seed)
        method = velo_method+'_'+gene_set_name
        split1 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
        split2 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
        df['split'+str(split_seed)] = compute_cosine_similarity_union(split1,split2,method)[0]
        print(split_seed)
    return df

df = create_cos_sim_df_across_seed(velo_method='sct', split_seeds=[317,320,323,326,329], grid_seeds=[227,230,233,236,239])
df.to_csv(data_folder+'cos_sim_across_seeds_'+gene_set_prefix+'_sct.csv')
"""
print('################# all done')

