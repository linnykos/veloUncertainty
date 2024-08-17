import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'sct'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 9815 × 2000
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split2_v2.h5ad') # 

### common genes being filtered out
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 1237
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(split1.var.index),np.array(total.var.index)).shape[0])) # 1531
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(split2.var.index),np.array(total.var.index)).shape[0])) # 1519

split1.layers['spliced_original'] = split1.layers['spliced']
split1.layers['unspliced_original'] = split1.layers['unspliced']
split2.layers['spliced_original'] = split2.layers['spliced']
split2.layers['unspliced_original'] = split2.layers['unspliced']

plot_gene_correlation_between_splits(split1,split2,fig_folder=fig_folder,fig_path="split_cor/"+dataset_short+"_"+method+"_corr_between_splits_colorSame.png")

exit()
##############################################
####################### old version
import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime
from sklearn.metrics.pairwise import cosine_similarity

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/" 

split1 = sc.read_h5ad(data_folder+'v2_erythroid/sct/adata_ery_sct_seed317_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_erythroid/sct/adata_ery_sct_seed317_split2_v2.h5ad')
total = sc.read_h5ad(data_folder+'v2_erythroid/sct/adata_ery_sct_total_v2.h5ad')

### common genes being selected and used for velocity computation
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index))
common_genes.shape # 1304->1237 common genes
np.intersect1d(np.array(split1.var.index), np.array(total.var.index)).shape # 229 -> 1531
np.intersect1d(np.array(split2.var.index), np.array(total.var.index)).shape # 244 -> 1519

# get the counts of commonly selected genes
gene_names_split1 = split1.var.index.copy()
positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
positions_split1 = [positions_dict_split1[gene] for gene in common_genes]

gene_names_split2 = split2.var.index.copy()
positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
positions_split2 = [positions_dict_split2[gene] for gene in common_genes]

def compute_gene_correlation_between_splits(mat1,mat2):
    m1 = np.transpose(mat1.todense())
    m2 = np.transpose(mat2.todense())
    cor = []
    for i in range(m1.shape[0]):
        if i%1000==0:
            print(i)
        if np.sum(m1[i,:])==0 and np.sum(m2[i,:])==0:
            cor.append(np.nan)
        else:
            cor.append(np.corrcoef(m1[i,:],m2[i,:])[0,1])
    cor = np.array(cor)
    return cor

# split1.layers['spliced'][:,positions_split1] is cell-by-gene, (9815, 1237)
cor_sct = compute_gene_correlation_between_splits(split1.layers['spliced'][:,positions_split1],split2.layers['spliced'][:,positions_split2])
cor_sct[~np.isnan(cor_sct)].shape
np.quantile(cor_sct[~np.isnan(cor_sct)],[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

## plot correlation between splits (all in one color)
plt.clf()
plt.scatter(range(len(cor_sct[~np.isnan(cor_sct)])), cor_sct[~np.isnan(cor_sct)],color='royalblue',alpha=0.5)
plt.title("Correlation of gene expr between splits (split_seed=317),Ngenes="+str(common_genes.shape[0]))
plt.xlabel("genes")
plt.ylabel("correlation")
#plt.savefig(fig_folder+"sct/split_cor/corr_between_splits_colorSame.png") 
plt.clf()

