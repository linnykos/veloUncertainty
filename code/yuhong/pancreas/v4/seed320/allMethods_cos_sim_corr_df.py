import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr, spearmanr

def compute_cosine_similarity_union(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    velo_split1 = pd.DataFrame(adata_split1.layers['velocity'], columns=velo_genes_split1)
    velo_split2 = pd.DataFrame(adata_split2.layers['velocity'], columns=velo_genes_split2)
    if method=='scv':
        velo_genes_split1 = velo_genes_split1[~np.isnan(velo_split1.loc[0])] #adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
        velo_genes_split2 = velo_genes_split2[~np.isnan(velo_split2.loc[0])] #adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    union_genes_velo = np.union1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Size of the union of genes for velocity computation in splits = '+str(union_genes_velo.shape[0])) 
    Nrow = adata_split1.shape[0]
    velo_df1 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split1:
        velo_df1[gene] = velo_split1[gene]
    velo_df2 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split2:
        velo_df2[gene] = velo_split2[gene]
    cos_sim = np.diag(cosine_similarity(velo_df1,velo_df2))
    return cos_sim, union_genes_velo.shape[0]

def compute_cosine_similarity_intersect(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    if method=="scv":
        velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
        velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0])) 
    velo_df1 = pd.DataFrame(adata_split1.layers['velocity'], columns=adata_split1.var.index.tolist())
    velo_df2 = pd.DataFrame(adata_split2.layers['velocity'], columns=adata_split2.var.index.tolist())
    cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))
    return cos_sim, common_genes_velocity.shape[0] # return cosine similarity and number of common genes in velocity computation

def compute_df_cos_sim_corr_across_methods(method, dataset_long, dataset_short, type, data_folder, seeds=[317,320,323,326,329]):
    outputAdded = ''
    if ((method=='sct') | ('velovi' in method)): outputAdded = '_outputAdded'
    for i in range(len(seeds)):
        seed = seeds[i]
        s1 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
        s2 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')
        if ('u' in type): 
            cos_sim = compute_cosine_similarity_union(s1,s2,method)[0]
            df[dataset_short+'_'+method+'_'+str(seed)] = cos_sim
        elif ('i' in type): 
            cos_sim = compute_cosine_similarity_intersect(s1,s2,method)[0]
            df[dataset_short+'_'+method+'_'+str(seed)] = cos_sim

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
dataset_long='pancreas'
dataset_short='pan'
df = pd.DataFrame()
for method in ['scv','utv','sct','velovi','velovi_woprep']:
    print(method)
    compute_df_cos_sim_corr_across_methods(method=method, dataset_long=dataset_long, dataset_short=dataset_short, 
                                           type='u', data_folder=data_folder, seeds=[317,320,323,326,329])

df.to_csv(data_folder+'v4_'+dataset_long+'/'+dataset_short+'_cos_sim_corr.csv')
