dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

total = ad.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 

# seed317
split_seed = 317
split1 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
split2 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')

genes_GPC = pd.read_csv(data_folder+'v4_'+dataset_long+'/glf_GPC.csv')
genes_keep = np.intersect1d(total.var.index, genes_GPC)

total_GPC = total[:, total.var.index.isin(genes_keep)]

split1_GPC = split1[:, split1.var.index.isin(genes_keep)]
split2_GPC = split2[:, split2.var.index.isin(genes_keep)]

total_GPC.write(data_folder+'v4_greenleaf/glf_total_GPC.h5ad')
split1_GPC.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_GPC.h5ad')
split2_GPC.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_GPC.h5ad')

# more seeds
total = ad.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 
genes_GPC = pd.read_csv(data_folder+'v4_'+dataset_long+'/glf_GPC.csv')
genes_keep = np.intersect1d(total.var.index, genes_GPC) # 184 genes

for split_seed in [320,323,326,329]:
    split1 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes.h5ad')
    split2 = sc.read(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes.h5ad')
    split1_GPC = split1[:, split1.var.index.isin(genes_keep)]
    split2_GPC = split2[:, split2.var.index.isin(genes_keep)]
    split1_GPC.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split1_GPC.h5ad')
    split2_GPC.write(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'_'+dataset_short+'_split2_GPC.h5ad')
    print(str(split_seed)+' done!')
