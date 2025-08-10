gene_set_name = 'GPC'
dataset_short = 'glf'
dataset_long = 'greenleaf'

import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
df_cos_sim = pd.DataFrame()
method_prefix = 'scv'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]

method_prefix = 'utv'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_GPC_utv.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_GPC_utv.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]

method_prefix = 'sct'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4_outputAdded.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]

method_prefix = 'velovi'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_GPC_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_GPC_outputAdded.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]

method_prefix = 'velovi_woprep'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_GPC_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_GPC_outputAdded.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]

df_cos_sim.to_csv(data_folder+'cos_sim_5seeds_glf_GPC.csv')

#########################################
dataset_short = 'glf'
dataset_long = 'greenleaf'

import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
df_cos_sim = pd.DataFrame()
gene_set_name_prefix = 'nGPCgrid'

method_prefix = 'scv'
for split_seed in [317,320,323,326,329]:
    grid_seed = split_seed-90
    gene_set_name = gene_set_name_prefix+str(grid_seed)
    method = method_prefix + '_' + gene_set_name
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]

method_prefix = 'utv'
for split_seed in [317,320,323,326,329]:
    grid_seed = split_seed-90
    gene_set_name = gene_set_name_prefix+str(grid_seed)
    method = method_prefix + '_' + gene_set_name
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]


# sct
dataset_short = 'glf'
dataset_long = 'greenleaf'
import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_sct import *

def compute_sct_avg_velocity(tnode,timesteps):
    v_shape = tnode.adata.shape
    v = np.zeros(v_shape)
    for t in timesteps:
        v += compute_sctour_velocity(tnode, timestep=t)
    return v/len(timesteps)

def sct_add_output(adata, tnode, timesteps):
    adata.layers['velocity'] = compute_sct_avg_velocity(tnode, timesteps)
    get_umap_sct(adata)

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
gene_set_name_prefix = 'nGPCgrid'
df_cos_sim = pd.read_csv(data_folder+'cos_sim_5seeds_glf_nGPCgrid.csv')

import sct
method_prefix = 'sct'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
    tnode_prefix = 'tnode_'+dataset_short+'_'+method
    tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
    tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')
    timesteps=[i/50 for i in range(1,11)]
    sct_add_output(split1, tnode_split1, timesteps)
    sct_add_output(split2, tnode_split2, timesteps)
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]



method_prefix = 'velovi'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_GPC_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_GPC_outputAdded.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]

method_prefix = 'velovi_woprep'
method = method_prefix + '_' + gene_set_name
for split_seed in [317,320,323,326,329]:
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_GPC_outputAdded.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_GPC_outputAdded.h5ad')
    df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]


df_cos_sim.to_csv(data_folder+'cos_sim_5seeds_glf_nGPCgrid.csv')


######################################
import os

os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
