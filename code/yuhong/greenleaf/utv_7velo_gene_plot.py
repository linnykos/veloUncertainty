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


###################################
dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'utv'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)

plot_mean_by_var(total.layers['velocity'], path=fig_folder+'gene_mean_var_total.png')
plot_mean_by_var(split1.layers['velocity'], path=fig_folder+'gene_mean_var_split1.png')
plot_mean_by_var(split2.layers['velocity'], path=fig_folder+'gene_mean_var_split2.png')

idx_maxvar = np.where(np.var(total.layers['velocity'], axis=0)>1000)[0][0]
plot_mean_by_var(np.delete(total.layers['velocity'], idx_maxvar, axis=1), 
                 path=fig_folder+'gene_mean_var_total_excludeMaxVar.png')


###################################
dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'utv_GPC'
split_seed=317

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/glf_total_GPC_utv.h5ad')
split1 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/seed'+str(split_seed)+'_glf_split1_GPC_utv.h5ad')
split2 = sc.read(data_folder+'seed'+str(split_seed)+'/'+method+'/seed'+str(split_seed)+'_glf_split2_GPC_utv.h5ad')

plot_mean_by_var(total.layers['velocity'], path=fig_folder+'gene_mean_var_total.png')
plot_mean_by_var(split1.layers['velocity'], path=fig_folder+'gene_mean_var_split1.png')
plot_mean_by_var(split2.layers['velocity'], path=fig_folder+'gene_mean_var_split2.png')

idx_maxvar = np.where(np.var(total.layers['velocity'], axis=0)>1000)[0][0]
plot_mean_by_var(np.delete(total.layers['velocity'], idx_maxvar, axis=1), 
                 path=fig_folder+'gene_mean_var_total_excludeMaxVar.png')

