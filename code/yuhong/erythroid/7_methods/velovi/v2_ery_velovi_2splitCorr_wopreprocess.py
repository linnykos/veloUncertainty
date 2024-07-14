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
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+'/wopreprocess/'

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/wopreprocess/adata_'+dataset_short+'_'+method+'_total_wopreprocess_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/wopreprocess/adata_'+dataset_short+'_'+method+'_split1_wopreprocess_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/wopreprocess/adata_'+dataset_short+'_'+method+'_split2_wopreprocess_v2.h5ad') # 

common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) # 1481
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 1481
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(split1.var.index),np.array(total.var.index)).shape[0])) # 1335
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(split2.var.index),np.array(total.var.index)).shape[0])) # 1328

plot_gene_correlation_between_splits(split1,split2,fig_path="split_cor/corr_between_splits_wopreprocess_colorSame.png",fig_folder=fig_folder)
"""
Spliced:
Number of valid values = 1481
Quantiles: [-1.19496360e-01 -1.01188167e-02 -6.07830647e-03 -3.05135275e-03
 -1.91648490e-04  2.12453379e-03  4.44817648e-03  7.44277937e-03
  1.07204593e-02  1.62678028e-02  9.68137230e-01]
Unspliced:
Number of valid values = 1481
Quantiles: [-0.02506132 -0.0098991  -0.00711452 -0.00508334 -0.00323733 -0.00147453
  0.00146106  0.00441732  0.00844299  0.01487345  0.50428819]
"""
