import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'velovi'
dataset_short = 'pan'
dataset_long = 'pancreas'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split1_v2.h5ad") # 3696 × 630
split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split2_v2.h5ad") # 3696 × 635
total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_total_v2.h5ad") # 3696 × 1074

### common genes being filtered out
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 521
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(split1.var.index),np.array(total.var.index)).shape[0])) # 584
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(split2.var.index),np.array(total.var.index)).shape[0])) # 583


plot_gene_correlation_between_splits(split1,split2,fig_folder=fig_folder,fig_path="/split_cor/"+dataset_short+"_"+method+"_corr_between_splits_colorSame.png")

"""
Spliced:
Number of valid values = 521
Quantiles: [-0.09006871 -0.0207196  -0.01447961 -0.01095189 -0.00638946 -0.00254371
  0.00156846  0.00575352  0.01151605  0.02082381  0.99589182]
Unspliced:
Number of valid values = 521
Quantiles: [-0.17162174 -0.02038644 -0.01528691 -0.01070023 -0.00644284 -0.0035595
  0.00099859  0.00658715  0.01401507  0.030669    0.84430089]
"""
#### wopreprocess
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'velovi'
dataset_short = 'pan'
dataset_long = 'pancreas'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+'/'+method+"/wopreprocess"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_v2.h5ad") # 3696 × 630
split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_v2.h5ad") # 3696 × 635
total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_v2.h5ad") # 3696 × 1074

### common genes being filtered out
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index)) 
print('Number of overlapped genes being filtered out in 2 splits = '+str(common_genes.shape[0])) # 1508
print('Number of overlapped genes being filtered out in split1 and total = '+str(np.intersect1d(np.array(split1.var.index),np.array(total.var.index)).shape[0])) # 1292
print('Number of overlapped genes being filtered out in splitw and total = '+str(np.intersect1d(np.array(split2.var.index),np.array(total.var.index)).shape[0])) # 1282

plot_gene_correlation_between_splits(split1,split2,fig_folder=fig_folder,fig_path="/split_cor/"+dataset_short+"_"+method+"_corr_between_splits_wopreprocess_colorSame.png")
"""
Spliced:
Number of valid values = 1508
Quantiles: [-0.15675754 -0.02006499 -0.01349503 -0.0093632  -0.00509138 -0.00112424
  0.00314461  0.0076622   0.0138886   0.02340991  0.99589182]
Unspliced:
Number of valid values = 1508
Quantiles: [-0.17162174 -0.01805836 -0.01222145 -0.0084482  -0.00523765 -0.00235818
  0.00175201  0.00658988  0.01218917  0.02115083  0.84430089]
"""