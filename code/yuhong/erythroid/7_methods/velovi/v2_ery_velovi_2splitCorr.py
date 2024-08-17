import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'velovi'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/" 

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 9815 × 310
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 9815 × 119
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 9815 × 128

### common genes being filtered out
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) # 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 91
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(split1.var.index),np.array(total.var.index)).shape[0])) # 114
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(split2.var.index),np.array(total.var.index)).shape[0])) # 125

plot_gene_correlation_between_splits(split1,split2,fig_path=method+"/split_cor/corr_between_splits_colorSame.png",fig_folder=fig_folder)


# color by sparsity

# color by housekeeping genes




### common genes being selected and used for velocity computation
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index))
common_genes.shape #  common genes
np.intersect1d(np.array(split1.var.index), np.array(total.var.index)).shape # 
np.intersect1d(np.array(split2.var.index), np.array(total.var.index)).shape # 

# get the counts of commonly selected genes
gene_names_split1 = split1.var.index.copy()
positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
positions_split1 = [positions_dict_split1[gene] for gene in common_genes]

gene_names_split2 = split2.var.index.copy()
positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
positions_split2 = [positions_dict_split2[gene] for gene in common_genes]
