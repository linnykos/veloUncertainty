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

s1 = sc.read_h5ad(data_folder+'v2_pancreas/seed317_split1_allgenes.h5ad')
s2 = sc.read_h5ad(data_folder+'v2_pancreas/seed317_split2_allgenes.h5ad')

np.corrcoef(np.array(s1.layers['spliced'].todense()[:,0].reshape(1,3696)),
            np.array(s2.layers['spliced'].todense()[:,0].reshape(1,3696)))
np.transpose(s1.layers['spliced'])

np.corrcoef(np.transpose(np.array(s1.layers['spliced'].todense()[:,0])),
            np.transpose(np.array(s2.layers['spliced'].todense()[0,:])))

np.corrcoef(np.transpose(s1.layers['spliced'].todense())[0,:],np.transpose(s2.layers['spliced'].todense())[0,:])

cor = [np.corrcoef(np.transpose(s1.layers['spliced'].todense())[i,:],np.transpose(s2.layers['spliced'].todense())[i,:]) for i in range(s1.layers['spliced'].shape[1])]

cor = []
S1 = np.transpose(s1.layers['spliced'].todense())
S2 = np.transpose(s2.layers['spliced'].todense())
for i in range(S1.shape[0]):
    if i%100==0:
        print(i)
    if np.sum(S1[i,:])==0 and np.sum(S2[i,:])==0:
        cor.append(np.nan)
    else:
        cor.append(np.corrcoef(S1[i,:],S2[i,:])[0,1])
cor = np.array(cor)
cor[~np.isnan(cor)].shape
np.quantile(cor[~np.isnan(cor)],[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])


split1 = sc.read_h5ad(data_folder+'v2_pancreas/sct/adata_pan_sct_seed317_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_pancreas/sct/adata_pan_sct_seed317_split2_v2.h5ad')
total = sc.read_h5ad(data_folder+'v2_pancreas/sct/adata_pan_sct_total_v2.h5ad')

common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) # 1335 common genes
np.intersect1d(np.array(split1.var.index), np.array(total.var.index)).shape # 1585
np.intersect1d(np.array(split2.var.index), np.array(total.var.index)).shape # 1571

velo_df1 = pd.DataFrame(split1.layers['velocity'], columns=split1.var.index.tolist())
velo_df2 = pd.DataFrame(split2.layers['velocity'], columns=split2.var.index.tolist())
cos_sim = np.diag(cosine_similarity(velo_df1[common_genes],velo_df2[common_genes]))




