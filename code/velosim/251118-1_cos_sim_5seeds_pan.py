import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

dataset_long = 'pancreas'
dataset_short = 'pan'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

method = 'scv'
df = pd.DataFrame(columns=[317,320,323,326,329])
for split_seed in [317,320,323,326,329]:
    adata1 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
    adata2 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
    df[split_seed] = compute_cosine_similarity_union(adata_split1=adata1,adata_split2=adata2,method=method)[0]

df.to_csv(data_folder+'v4_'+dataset_long+'/cos_sim_5seeds_'+dataset_short+'_'+method+'.csv')

method = 'utv'
df = pd.DataFrame(columns=[317,320,323,326,329])
for split_seed in [317,320,323,326,329]:
    adata1 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
    adata2 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
    df[split_seed] = compute_cosine_similarity_union(adata_split1=adata1,adata_split2=adata2,method=method)[0]

df.to_csv(data_folder+'v4_'+dataset_long+'/cos_sim_5seeds_'+dataset_short+'_'+method+'.csv')

