dataset_long = 'greenleaf'
dataset_short = 'glf'
split_seed=317

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *
from v4_functions import *


def compute_velocity_union_df(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    velo_split1 = pd.DataFrame(adata_split1.layers['velocity'], columns=velo_genes_split1)
    velo_split2 = pd.DataFrame(adata_split2.layers['velocity'], columns=velo_genes_split2)
    if 'scv' in method:
        velo_genes_split1 = velo_genes_split1[~np.isnan(velo_split1.loc[0])] 
        velo_genes_split2 = velo_genes_split2[~np.isnan(velo_split2.loc[0])] 
    union_genes_velo = np.union1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Size of the union of genes for velocity computation in splits = '+str(union_genes_velo.shape[0])) 
    Nrow = adata_split1.shape[0]
    velo_df1 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split1:
        velo_df1[gene] = velo_split1[gene]
    velo_df2 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split2:
        velo_df2[gene] = velo_split2[gene]
    return velo_df1, velo_df2

# scv
method = 'scv'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split1_v4.h5ad')
split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split2_v4.h5ad')
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)
np.mean( [cos_sim[i,i]>np.quantile(cos_sim[i], .95) for i in range(df1.shape[0])] ) # proportion lying on the 5% tail

method = 'scv_GPC'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
split1 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
split2 = sc.read_h5ad(data_folder+'adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)

# utv
method = 'utv'
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)
np.mean( [cos_sim[i,i]>np.quantile(cos_sim[i], .95) for i in range(df1.shape[0])] ) # proportion lying on the 5% tail

method = 'utv_GPC'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
split1 = sc.read_h5ad(data_folder+'seed317_glf_split1_GPC_utv.h5ad')
split2 = sc.read_h5ad(data_folder+'seed317_glf_split2_GPC_utv.h5ad')
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)

# sct
method = 'sct'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
df1, df2 = compute_velocity_union_df(split1,split2,method='sct')
cos_sim = cosine_similarity(df1, df2)
np.mean( [cos_sim[i,i]>np.quantile(cos_sim[i], .95) for i in range(df1.shape[0])] ) # proportion lying on the 5% tail

method = 'sct_GPC'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
df1, df2 = compute_velocity_union_df(split1,split2,method='sct')
cos_sim = cosine_similarity(df1, df2)

# velovi
method='velovi'
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)
np.mean( [cos_sim[i,i]>np.quantile(cos_sim[i], .95) for i in range(df1.shape[0])] ) # proportion lying on the 5% tail

method='velovi_GPC'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
split1 = sc.read(data_folder+'adata_glf_velovi_GPC_split1_GPC_outputAdded.h5ad')
split2 = sc.read(data_folder+'adata_glf_velovi_GPC_split2_GPC_outputAdded.h5ad')
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)

# velovi_woprep
method='velovi_woprep'
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)
np.mean( [cos_sim[i,i]>np.quantile(cos_sim[i], .95) for i in range(df1.shape[0])] ) # proportion lying on the 5% tail

method='velovi_woprep_GPC'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
split1 = sc.read(data_folder+'adata_glf_'+method+'_split1_GPC_outputAdded.h5ad')
split2 = sc.read(data_folder+'adata_glf_'+method+'_split2_GPC_outputAdded.h5ad')
df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)


from sklearn.metrics.pairwise import cosine_similarity
for celltype_cat in total.obs.cluster_name.cat.categories:
    # do the shuffling per cell within a celltype
    cell_idx = np.where(total.obs.cluster_name==celltype_cat)[0]
    shuf_cos_sim_cells = cosine_similarity(df1.iloc[cell_idx], df2.iloc[cell_idx])
    shuf_quantiles = np.array([np.mean(shuf_cos_sim_cells[i,:]>=shuf_cos_sim_cells[i,i]) for i in range(shuf_cos_sim_cells.shape[0])])
    print(celltype_cat + ' '+ str( np.round(np.mean(shuf_quantiles<=.05),4) ) )

df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)
off_diagonal = cos_sim[~np.eye(cos_sim.shape[0], dtype=bool)]
np.quantile(off_diagonal, [0.,.25,.5,.75,1.])
offdiag_q95 = np.quantile(off_diagonal, .95) 
np.diag(cos_sim)
np.sum(np.diag(cos_sim)>offdiag_q95)
np.mean(np.diag(cos_sim)>offdiag_q95) 

df1, df2 = compute_velocity_union_df(split1,split2,method=method)
cos_sim = cosine_similarity(df1, df2)
np.mean( [cos_sim[i,i]>np.quantile(cos_sim[i], .95) for i in range(df1.shape[0])] ) # proportion lying on the 5% tail
# scv: 0.139 -> 0.064
# utv: 0.1406 -> 0.135
# sct: 0.1038 -> 0.2658
# velovi: 0.5004 -> 0.1339
# velovi_woprep: 0.2388 -> 0.1275

csq95 = np.quantile(cos_sim, .95)
csq95
np.mean( [cos_sim[i,i]>csq95 for i in range(df1.shape[0])] )
np.mean( [ np.diag(cos_sim)>csq95 ] )
# scv: 0.1369 (q95=0.6292) -> 0.0618 (q95=0.9458)
# utv: 0.0816 (q95=0.9425) -> 0.1636 (q95=0.9399)
# sct: 0.0271 (q95=0.8568) -> 0.1587 (q95=0.9811)
# velovi: 0.2716 (q95=0.6249) -> 0.1268 (q95=0.6550)
# velovi_woprep: 0.1561 (q95=0.8934) -> 0.1216 (q95=0.9471)

