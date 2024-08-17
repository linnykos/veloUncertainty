import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = "utv"
dataset_short = "pan"
dataset_long = "pancreas"

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/" 

# cd ~/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup16_v2_utv
total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

### common genes being filtered out
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 1508
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(split1.var.index),np.array(total.var.index)).shape[0])) # 1292
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(split2.var.index),np.array(total.var.index)).shape[0])) # 1282


plot_gene_correlation_between_splits(split1,split2,fig_folder=fig_folder,fig_path=method+"/split_cor/"+dataset_short+"_"+method+"_corr_between_splits_colorSame.png")

