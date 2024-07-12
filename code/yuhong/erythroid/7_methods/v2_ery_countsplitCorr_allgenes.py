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

# correlation between splits
split1_allgenes = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split1_allgenes.h5ad')
split2_allgenes = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split2_allgenes.h5ad')

cor = []
S1 = np.transpose(split1_allgenes.layers['spliced'].todense())
S2 = np.transpose(split2_allgenes.layers['spliced'].todense())
for i in range(S1.shape[0]):
    if i%1000==0:
        print(i)
    if np.sum(S1[i,:])==0 and np.sum(S2[i,:])==0:
        cor.append(np.nan)
    else:
        cor.append(np.corrcoef(S1[i,:],S2[i,:])[0,1])
cor = np.array(cor)
cor[~np.isnan(cor)].shape
np.quantile(cor[~np.isnan(cor)],[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

## plot correlation between splits (all in one color)
plt.clf()
plt.scatter(range(len(cor[~np.isnan(cor)])), cor[~np.isnan(cor)],color='royalblue',alpha=0.5)
plt.title("Correlation of gene expr between splits (split_seed=317)")
plt.xlabel("genes (with nonzero expr)")
plt.ylabel("correlation")
plt.savefig(fig_folder+"corr_between_splits_allgenes.png") 
plt.clf()
