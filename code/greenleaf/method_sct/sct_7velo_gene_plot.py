import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

import numpy as np
import matplotlib.pyplot as plt
## mean~variance for each gene
def plot_mean_by_var(velo, path):
    means = np.mean(velo, axis=0)
    variances = np.var(velo, axis=0)
    plt.scatter(means, variances)
    plt.xlabel('Mean')
    plt.ylabel('Variance')
    plt.title('Mean vs variance per gene')
    plt.grid(True)
    plt.savefig(path)
    plt.clf()

#####################################
split_seed = 317
method = 'sct'
dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"

adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

plot_mean_by_var(total.layers['velocity'], path=fig_folder+'gene_mean_var_total.png')
plot_mean_by_var(split1.layers['velocity'], path=fig_folder+'gene_mean_var_split1.png')
plot_mean_by_var(split2.layers['velocity'], path=fig_folder+'gene_mean_var_split2.png')

#####################################
split_seed = 317
method = 'sct_GPC'
dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"

adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

plot_mean_by_var(total.layers['velocity'], path=fig_folder+'gene_mean_var_total.png')
plot_mean_by_var(split1.layers['velocity'], path=fig_folder+'gene_mean_var_split1.png')
plot_mean_by_var(split2.layers['velocity'], path=fig_folder+'gene_mean_var_split2.png')

