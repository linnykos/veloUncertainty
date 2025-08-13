import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
import os
import re

dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

## nGPC
gene_set_name_prefix = 'nGPCblk'
df_cos_sim = pd.DataFrame()

for method_prefix in ['scv','utv','sct','velovi','velovi_woprep']:
    for i in range(5):
        split_seed = [317,320,323,326,329][i]
        grid_seed = [227,230,233,236,239][i]
        gene_set_name = gene_set_name_prefix + str(grid_seed)
        method = method_prefix + '_' + gene_set_name
        adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
        if 'scv' in method or 'utv' in method:
            path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1', s)][0]
            path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2', s)][0]
        else: 
            path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1.*outputAdded', s)][0]
            path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2.*outputAdded', s)][0]
        split1 = sc.read_h5ad( path_split1 )
        split2 = sc.read_h5ad( path_split2 )
        print(path_split1)
        print(path_split2)
        df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]
    print(method_prefix+' done.')

df_cos_sim['cell'] = split1.obs.index
df_cos_sim['celltype'] = np.array(split1.obs.cluster_name)
df_cos_sim.to_csv(data_folder+'cos_sim_5seeds_'+dataset_short+'_'+gene_set_name_prefix+'.csv')

# df_cos_sim = pd.read_csv(data_folder+'cos_sim_5seeds_'+dataset_short+'_'+gene_set_name_prefix+'.csv')

## Mark
gene_set_name = 'GPC'
df_cos_sim = pd.DataFrame()

for method_prefix in ['scv','utv','sct','velovi','velovi_woprep']:
    for i in range(5):
        split_seed = [317,320,323,326,329][i]
        method = method_prefix + '_' + gene_set_name
        adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
        if 'scv' in method or 'utv' in method:
            path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1', s)][0]
            path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2', s)][0]
        else: 
            path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1.*outputAdded', s)][0]
            path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2.*outputAdded', s)][0]
        split1 = sc.read_h5ad( path_split1 )
        split2 = sc.read_h5ad( path_split2 )
        print(path_split1)
        print(path_split2)
        df_cos_sim[method_prefix+'_seed'+str(split_seed)] = compute_cosine_similarity_union(adata_split1=split1, adata_split2=split2, method=method_prefix)[0]
    print(method_prefix+' done.')

df_cos_sim['cell'] = split1.obs.index
df_cos_sim['celltype'] = np.array(split1.obs.cluster_name)

df_cos_sim.to_csv(data_folder+'cos_sim_5seeds_'+dataset_short+'_'+gene_set_name+'.csv')

