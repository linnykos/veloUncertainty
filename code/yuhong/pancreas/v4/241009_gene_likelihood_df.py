import pandas as pd
import numpy as np
import scanpy as sc

method = 'scv'

s1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed317/'+method+'/adata_pan_'+method+'_split1_v4.h5ad')
s2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed317/'+method+'/adata_pan_'+method+'_split2_v4.h5ad')
total = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed317/scv/adata_pan_scv_total_v4.h5ad')
s1.var['fit_likelihood']

genes = total.var.index
for seed in [317,320,323,326,329]:
    s1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed'+str(seed)+'/scv/adata_pan_scv_split1_v4.h5ad')
    s2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed'+str(seed)+'/scv/adata_pan_scv_split2_v4.h5ad')
    genes = np.union1d(genes, np.union1d(s1.var.index, s2.var.index))

genes = np.sort(genes)

df_lik = pd.DataFrame()
df_lik['gene_name'] = genes
df_lik['total_highly_variable'] = np.nan
df_lik['total_velocity_genes'] = np.nan
df_lik['total_lik'] = np.nan
for gene in total.var.index:
    idx = np.where(df_lik['gene_name'] == gene)[0][0]
    df_lik.loc[idx, 'total_highly_variable'] = total.var['highly_variable'][gene]
    df_lik.loc[idx, 'total_velocity_genes'] = total.var['velocity_genes'][gene]
    df_lik.loc[idx, 'total_lik'] = total.var['fit_likelihood'][gene]
    

# np.sum(df_lik['total_highly_variable']==1)
# np.sum(total.var['highly_variable']==False)
# np.sum(total.var['velocity_genes']==False)
# np.sum(df_lik['total_velocity_genes']==1)

for seed in [317,320,323,326,329]:
    s1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed'+str(seed)+'/scv/adata_pan_scv_split1_v4.h5ad')
    s2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed'+str(seed)+'/scv/adata_pan_scv_split2_v4.h5ad')
    df_lik[str(seed)+'s1_highly_variable'] = np.nan
    df_lik[str(seed)+'s1_velocity_genes'] = np.nan
    df_lik[str(seed)+'s1_lik'] = np.nan
    df_lik[str(seed)+'s2_highly_variable'] = np.nan
    df_lik[str(seed)+'s2_velocity_genes'] = np.nan
    df_lik[str(seed)+'s2_lik'] = np.nan
    for gene in s1.var.index:
        idx = np.where(df_lik['gene_name'] == gene)[0][0]
        df_lik.loc[idx, str(seed)+'s1_lik'] = s1.var['fit_likelihood'][gene]
        df_lik.loc[idx, str(seed)+'s1_highly_variable'] = s1.var['highly_variable'][gene]
        df_lik.loc[idx, str(seed)+'s1_velocity_genes'] = s1.var['velocity_genes'][gene]
    for gene in s2.var.index:
        idx = np.where(df_lik['gene_name'] == gene)[0][0]
        df_lik.loc[idx, str(seed)+'s2_lik'] = s2.var['fit_likelihood'][gene]
        df_lik.loc[idx, str(seed)+'s2_highly_variable'] = s2.var['highly_variable'][gene]
        df_lik.loc[idx, str(seed)+'s2_velocity_genes'] = s2.var['velocity_genes'][gene]

np.sum(df_lik['317s1_highly_variable']==1)
np.sum(df_lik['317s1_highly_variable']==0)
np.sum(df_lik['317s1_velocity_genes']==0)
np.sum(np.isnan(df_lik['317s1_highly_variable']))
np.sum(~np.isnan(df_lik['317s1_highly_variable']))

# check
for seed in [317,320,323,326,329]:
    print(seed)
    print(np.sum(df_lik[str(seed)+'s1_highly_variable']==1)+np.sum(df_lik[str(seed)+'s1_highly_variable']==0),
          np.sum(df_lik[str(seed)+'s1_highly_variable']==1),np.sum(df_lik[str(seed)+'s1_highly_variable']==0),
          np.sum(df_lik[str(seed)+'s1_velocity_genes']==1),np.sum(df_lik[str(seed)+'s1_velocity_genes']==0))

raw = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/endocrinogenesis_day15.h5ad')
allgenes = raw.var.index
overdisp_S = pd.read_csv('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/pan_overdisp_S.csv')
overdisp_U = pd.read_csv('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/pan_overdisp_U.csv')

for gene in genes:
    idx = np.where(df_lik['gene_name'] == gene)[0][0]
    idx_allgenes = np.where(allgenes==gene)[0][0]
    df_lik.loc[idx, 'overdisp_S'] = overdisp_S['x'][idx_allgenes]
    df_lik.loc[idx, 'overdisp_U'] = overdisp_U['x'][idx_allgenes]

for gene in genes:
    idx = np.where(df_lik['gene_name'] == gene)[0][0]
    idx_allgenes = np.where(allgenes==gene)[0][0]
    df_lik.loc[idx, 'nonzero_S'] = np.mean(raw.layers['spliced'][:,idx_allgenes] > 0)
    df_lik.loc[idx, 'nonzero_U'] = np.mean(raw.layers['unspliced'][:,idx_allgenes] > 0)

df_lik.to_csv('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/pancreas_gene_likelihood.csv')

df_lik[['317s1_highly_variable', '317s1_lik', '317s2_highly_variable',
       '317s2_lik', '320s1_highly_variable', '320s1_lik',
       '320s2_highly_variable', '320s2_lik', '323s1_highly_variable',
       '323s1_lik', '323s2_highly_variable', '323s2_lik',
       '326s1_highly_variable', '326s1_lik', '326s2_highly_variable',
       '326s2_lik', '329s1_highly_variable', '329s1_lik',
       '329s2_highly_variable', '329s2_lik']].iloc[np.where(df_lik['total_lik']>.4)[0]]

df_lik[['317s1_highly_variable', '317s1_lik', '317s2_highly_variable', '317s2_lik']].iloc[np.where(df_lik['total_lik']>.4)[0]]

df_lik['317s1_highly_variable'].iloc[np.where(df_lik['total_lik']>.4)[0]].values
df_lik['317s2_highly_variable'].iloc[np.where(df_lik['total_lik']>.4)[0]].values

for seed in [317,320,323,326,329]:
    print(seed)
    col_name = [str(seed)+'s1_highly_variable', str(seed)+'s2_highly_variable']
    print(np.sum(df_lik[col_name].iloc[np.where(df_lik['total_lik']>.4)[0]].values))
