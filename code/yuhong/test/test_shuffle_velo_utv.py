import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'utv'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

######################################################
dataset_long = 'erythroid'
dataset_short = 'ery'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

"""
np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
= (0.64563376, 0.095295645)

np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
= (0.042651482, 0.008892522)

(np.mean(cos_sim)-np.mean(v2s_mean))/np.sqrt(np.var(v2s_mean))
= 93.36633

(np.median(cos_sim)-np.mean(v2s_median))/np.sqrt(np.var(v2s_median))
= 3.701508
"""
######################################################
dataset_long = 'pancreas'
dataset_short = 'pan'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

"""
>>> np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
(0.39516446, 0.034671333)
>>> 
>>> np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
(0.2585881, 0.009260354)
>>> 
>>> (np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))
65.73474
>>> 
>>> (np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))
13.045775
"""

######################################################
dataset_long = 'pancreasINC'
dataset_short = 'panINC'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

"""
>>> np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
(0.3420931, 0.02539015)
>>> np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
(0.19365247, 0.01110906)
>>> (np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))
62.583477
>>> (np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))
21.232706
"""
######################################################
dataset_long = 'larry'
dataset_short = 'larry'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

"""
>>> np.round(np.mean(cos_sim),4), np.round(np.median(cos_sim),4)
(0.8058, 0.9813)

>>> np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
(0.81706834, 0.0038185308)
>>> np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
(0.32629508, 0.0026658222)
>>> (np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))
an)) / np.sqrt(np.var(v2s_median))
179.85513
>>> (np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))
43.00862

np.mean(v2s_mean),np.median(v2s_mean) # (0.3254157, 0.3251)
np.mean(v2s_median),np.median(v2s_median) # (0.81595486, 0.8161)
np.min(v2s_mean),np.max(v2s_mean) # (0.3216, 0.3305)
np.min(v2s_median),np.max(v2s_median) # (0.8111, 0.8236)
"""

v1 = split1.layers['velocity'].copy()
v2 = split2.layers['velocity'].copy()
v1 = pd.DataFrame(v1)
v2 = pd.DataFrame(v2)
v1.columns = split1.var.index
v2.columns = split2.var.index
common_genes_velocity = np.intersect1d(np.array(v1.columns), np.array(v2.columns))
cos_sim = np.diag(cosine_similarity(v1[common_genes_velocity],v2[common_genes_velocity]))
np.round(np.mean(cos_sim),4), np.round(np.median(cos_sim),4)

v2s_mean = []
v2s_median = []
np.random.seed(1514)
for i in range(101):
    if i % 10==0: print(i)
    v2s = v2.sample(frac=1)
    cos_sim_s = np.diag(cosine_similarity(v1[common_genes_velocity],v2s[common_genes_velocity]))
    v2s_mean.append(np.round(np.mean(cos_sim_s),4))
    v2s_median.append(np.round(np.median(cos_sim_s),4))

np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 

np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 

(np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))

(np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))


np.mean(v2s_mean),np.median(v2s_mean) 
np.mean(v2s_median),np.median(v2s_median) 
np.var(v2s_mean)
np.min(v2s_mean),np.max(v2s_mean)
np.min(v2s_median),np.max(v2s_median)

(np.mean(cos_sim)-np.mean(v2s_mean)) /np.sqrt(np.var(v2s_mean))
(np.median(cos_sim)-np.mean(v2s_median)) /np.sqrt(np.var(v2s_median))
"""
v1 = split1.layers['velocity'][:,~np.isnan(split1.layers['velocity'][0,:])]
v1_columns = split1.var.index[~np.isnan(split1.layers['velocity'][0,:])]
v2 = split2.layers['velocity'][:,~np.isnan(split2.layers['velocity'][0,:])]
v2_columns = split2.var.index[~np.isnan(split2.layers['velocity'][0,:])]
common_genes_velocity = np.intersect1d(np.array(v1_columns), np.array(v2_columns))

v1 = pd.DataFrame(v1)
v1.columns = v1_columns
v2 = pd.DataFrame(v2)
v2.columns = v2_columns
cos_sim = np.diag(cosine_similarity(v1[common_genes_velocity],v2[common_genes_velocity]))
np.round(np.mean(cos_sim),4), np.round(np.median(cos_sim),4)

np.random.seed(1226)
v2s = v2.sample(frac=1)
cos_sim_s = np.diag(cosine_similarity(v1[common_genes_velocity],v2s[common_genes_velocity]))
np.quantile(cos_sim_s, [0,.25,.5,.75,1.])
np.round(np.mean(cos_sim_s),4), np.round(np.median(cos_sim_s),4)
"""