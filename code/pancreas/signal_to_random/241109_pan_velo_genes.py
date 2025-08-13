import pandas as pd
import numpy as np
import scanpy as sc

method = 'utv'

#s1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed317/'+method+'/adata_pan_'+method+'_split1_v4.h5ad')
#s2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/seed317/'+method+'/adata_pan_'+method+'_split2_v4.h5ad')

# s1.var['fit_loss']
# s1.var['highly_variable']

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_pancreas/'

total = sc.read_h5ad(data_folder+'seed317/'+method+'/adata_pan_'+method+'_total_v4.h5ad')
genes = total.var.index
for seed in [317,320,323,326,329]:
    s1 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split1_v4.h5ad')
    s2 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split2_v4.h5ad')
    genes = np.union1d(genes, np.union1d(s1.var.index, s2.var.index))

genes = np.sort(genes)

df = pd.DataFrame()
df['gene_name'] = genes
df['total_highly_variable'] = np.nan
df['total_velocity_genes'] = np.nan
df['total_fit_loss'] = np.nan
df['total_fit_llf'] = np.nan
df['total_fit_sr2'] = np.nan
df['total_fit_ur2'] = np.nan
for gene in total.var.index:
    idx = np.where(df['gene_name'] == gene)[0][0]
    df.loc[idx, 'total_highly_variable'] = total.var['highly_variable'][gene]
    df.loc[idx, 'total_velocity_genes'] = total.var['velocity_genes'][gene]
    df.loc[idx, 'total_fit_loss'] = total.var['fit_loss'][gene]
    df.loc[idx, 'total_fit_llf'] = total.var['fit_llf'][gene]
    df.loc[idx, 'total_fit_sr2'] = total.var['fit_sr2'][gene]
    df.loc[idx, 'total_fit_ur2'] = total.var['fit_ur2'][gene]
    
for seed in [317,320,323,326,329]:
    # s1
    s1 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split1_v4.h5ad')
    df[str(seed)+'s1_highly_variable'] = np.nan
    df[str(seed)+'s1_velocity_genes'] = np.nan
    df[str(seed)+'s1_fit_loss'] = np.nan
    df[str(seed)+'s1_fit_llf'] = np.nan
    df[str(seed)+'s1_fit_sr2'] = np.nan
    df[str(seed)+'s1_fit_ur2'] = np.nan
    # s2
    s2 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split2_v4.h5ad')
    df[str(seed)+'s2_highly_variable'] = np.nan
    df[str(seed)+'s2_velocity_genes'] = np.nan
    df[str(seed)+'s2_fit_loss'] = np.nan
    df[str(seed)+'s2_fit_llf'] = np.nan
    df[str(seed)+'s2_fit_sr2'] = np.nan
    df[str(seed)+'s2_fit_ur2'] = np.nan
    for gene in s1.var.index:
        idx = np.where(df['gene_name'] == gene)[0][0]
        df.loc[idx, str(seed)+'s1_fit_loss'] = s1.var['fit_loss'][gene]
        df.loc[idx, str(seed)+'s1_fit_llf'] = s1.var['fit_llf'][gene]
        df.loc[idx, str(seed)+'s1_fit_sr2'] = s1.var['fit_sr2'][gene]
        df.loc[idx, str(seed)+'s1_fit_ur2'] = s1.var['fit_ur2'][gene]
        df.loc[idx, str(seed)+'s1_highly_variable'] = s1.var['highly_variable'][gene]
        df.loc[idx, str(seed)+'s1_velocity_genes'] = s1.var['velocity_genes'][gene]
    for gene in s2.var.index:
        idx = np.where(df['gene_name'] == gene)[0][0]
        df.loc[idx, str(seed)+'s2_fit_loss'] = s2.var['fit_loss'][gene]
        df.loc[idx, str(seed)+'s2_fit_llf'] = s2.var['fit_llf'][gene]
        df.loc[idx, str(seed)+'s2_fit_sr2'] = s2.var['fit_sr2'][gene]
        df.loc[idx, str(seed)+'s2_fit_ur2'] = s2.var['fit_ur2'][gene]
        df.loc[idx, str(seed)+'s2_highly_variable'] = s2.var['highly_variable'][gene]
        df.loc[idx, str(seed)+'s2_velocity_genes'] = s2.var['velocity_genes'][gene]


raw = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad')
allgenes = raw.var.index
overdisp_S = pd.read_csv(data_folder+'pan_overdisp_S.csv')
overdisp_U = pd.read_csv(data_folder+'pan_overdisp_U.csv')

for gene in genes:
    idx = np.where(df['gene_name'] == gene)[0][0]
    idx_allgenes = np.where(allgenes==gene)[0][0]
    df.loc[idx, 'overdisp_S'] = overdisp_S['x'][idx_allgenes]
    df.loc[idx, 'overdisp_U'] = overdisp_U['x'][idx_allgenes]

for gene in genes:
    idx = np.where(df['gene_name'] == gene)[0][0]
    idx_allgenes = np.where(allgenes==gene)[0][0]
    df.loc[idx, 'nonzero_S'] = np.mean(raw.layers['spliced'][:,idx_allgenes] > 0)
    df.loc[idx, 'nonzero_U'] = np.mean(raw.layers['unspliced'][:,idx_allgenes] > 0)

#df['317s1_highly_variable'].iloc[np.where(df['total_fit_loss']>.4)[0]].values
#df['317s2_highly_variable'].iloc[np.where(df['total_fit_loss']>.4)[0]].values


df.to_csv(data_folder+'pancreas_'+method+'_gene_fitloss.csv')


###################
import pandas as pd
import numpy as np
import scanpy as sc

method = 'velovi'
# s1.var['velocity_r2']
# s1.var['highly_variable']
# s1.var['velocity_genes']

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_pancreas/'

total = sc.read_h5ad(data_folder+'seed317/'+method+'/adata_pan_'+method+'_total_v4.h5ad')
genes = total.var.index
for seed in [317,320,323,326,329]:
    s1 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split1_v4.h5ad')
    s2 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split2_v4.h5ad')
    genes = np.union1d(genes, np.union1d(s1.var.index, s2.var.index))

genes = np.sort(genes)

df = pd.DataFrame()
df['gene_name'] = genes
df['total_highly_variable'] = np.nan
df['total_velocity_genes'] = np.nan
df['total_velocity_r2'] = np.nan
for gene in total.var.index:
    idx = np.where(df['gene_name'] == gene)[0][0]
    df.loc[idx, 'total_highly_variable'] = total.var['highly_variable'][gene]
    df.loc[idx, 'total_velocity_genes'] = total.var['velocity_genes'][gene]
    df.loc[idx, 'total_velocity_r2'] = total.var['velocity_r2'][gene]
   
for seed in [317,320,323,326,329]:
    s1 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split1_v4_outputAdded.h5ad')
    s2 = sc.read_h5ad(data_folder+'seed'+str(seed)+'/'+method+'/adata_pan_'+method+'_split2_v4_outputAdded.h5ad')
    df[str(seed)+'s1_highly_variable'] = np.nan
    df[str(seed)+'s1_velocity_genes'] = np.nan
    df[str(seed)+'s1_velocity_r2'] = np.nan
    df[str(seed)+'s2_highly_variable'] = np.nan
    df[str(seed)+'s2_velocity_genes'] = np.nan
    df[str(seed)+'s2_velocity_r2'] = np.nan
    for gene in s1.var.index:
        idx = np.where(df['gene_name'] == gene)[0][0]
        df.loc[idx, str(seed)+'s1_velocity_r2'] = s1.var['velocity_r2'][gene]
        df.loc[idx, str(seed)+'s1_highly_variable'] = s1.var['highly_variable'][gene]
        df.loc[idx, str(seed)+'s1_velocity_genes'] = s1.var['velocity_genes'][gene]
    for gene in s2.var.index:
        idx = np.where(df['gene_name'] == gene)[0][0]
        df.loc[idx, str(seed)+'s2_velocity_r2'] = s2.var['velocity_r2'][gene]
        df.loc[idx, str(seed)+'s2_highly_variable'] = s2.var['highly_variable'][gene]
        df.loc[idx, str(seed)+'s2_velocity_genes'] = s2.var['velocity_genes'][gene]


raw = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad')
allgenes = raw.var.index
overdisp_S = pd.read_csv(data_folder+'pan_overdisp_S.csv')
overdisp_U = pd.read_csv(data_folder+'pan_overdisp_U.csv')

for gene in genes:
    idx = np.where(df['gene_name'] == gene)[0][0]
    idx_allgenes = np.where(allgenes==gene)[0][0]
    df.loc[idx, 'overdisp_S'] = overdisp_S['x'][idx_allgenes]
    df.loc[idx, 'overdisp_U'] = overdisp_U['x'][idx_allgenes]

for gene in genes:
    idx = np.where(df['gene_name'] == gene)[0][0]
    idx_allgenes = np.where(allgenes==gene)[0][0]
    df.loc[idx, 'nonzero_S'] = np.mean(raw.layers['spliced'][:,idx_allgenes] > 0)
    df.loc[idx, 'nonzero_U'] = np.mean(raw.layers['unspliced'][:,idx_allgenes] > 0)


df.to_csv(data_folder+'pancreas_'+method+'_gene_velocityR2.csv')



