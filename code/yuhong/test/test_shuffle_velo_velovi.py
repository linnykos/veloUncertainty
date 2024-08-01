import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'velovi'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

def compute_cosine_similarity_shuffled(split1,split2,seed=1514,Niter=100):
    from sklearn.metrics.pairwise import cosine_similarity
    v1 = split1.layers['velocity'].copy()
    v2 = split2.layers['velocity'].copy()
    v1 = pd.DataFrame(v1)
    v2 = pd.DataFrame(v2)
    v1.columns = split1.var.index
    v2.columns = split2.var.index
    common_genes_velocity = np.intersect1d(np.array(v1.columns), np.array(v2.columns))
    cos_sim = np.diag(cosine_similarity(v1[common_genes_velocity],v2[common_genes_velocity]))
    print(np.round(np.mean(cos_sim),4), np.round(np.median(cos_sim),4))
    v2s_mean = []
    v2s_median = []
    np.random.seed(seed)
    for i in range(Niter+1):
        if i % 10==0: print(i)
        v2s = v2.sample(frac=1)
        cos_sim_s = np.diag(cosine_similarity(v1[common_genes_velocity],v2s[common_genes_velocity]))
        v2s_mean.append(np.round(np.mean(cos_sim_s),4))
        v2s_median.append(np.round(np.median(cos_sim_s),4))
    return v2s_mean,v2s_median

######################################################
dataset_long = 'erythroid'
dataset_short = 'ery'
celltype_label = 'celltype'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded_v2.h5ad')

v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,seed=1514)
"""
>>> np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
(0.39259705, 0.002341551)
>>> np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
(0.4190802, 0.0015593243)
>>> (np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))
77.1564
>>> (np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))
83.33109
"""
######################################################
dataset_long = 'pancreas'
dataset_short = 'pan'
celltype_label = 'clusters'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded_v2.h5ad')

v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,seed=1514)

"""
>>> np.round(np.mean(cos_sim),4), np.round(np.median(cos_sim),4)
(0.4501, 0.4531)
>>> np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
(0.23284158, 0.0024521668)
>>> np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
(0.26088315, 0.0018664864)
>>> (np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))
101.377365
>>> (np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))
89.8296
"""

######################################################
dataset_long = 'pancreasINC'
dataset_short = 'panINC'
celltype_label = 'clusters'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded_v2.h5ad')

v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2)

"""
>>> np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
(0.14694159, 0.0023807092)
>>> np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
(0.1834812, 0.0018501121)
>>> (np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))
144.11095
>>> (np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))
128.60753
"""
######################################################
dataset_long = 'larry'
dataset_short = 'larry'
celltype_label = 'state_info'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded_v2.h5ad')


v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2)
"""
np.mean(cos_sim)=0.58 
np.median(cos_sim)=0.5915
>>> np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 
(0.24432872, 0.0007919569)
>>> np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 
(0.2761871, 0.0005725846)

"""

np.mean(v2s_median), np.sqrt(np.var(v2s_median)) 

np.mean(v2s_mean), np.sqrt(np.var(v2s_mean)) 

#(np.mean(cos_sim)-np.mean(v2s_mean)) / np.sqrt(np.var(v2s_mean))
#(np.median(cos_sim)-np.mean(v2s_median)) / np.sqrt(np.var(v2s_median))
