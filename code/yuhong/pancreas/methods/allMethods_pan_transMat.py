import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_transMat import *

dataset_long = 'larry'
dataset_short = 'larry'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

## scv
method = 'scv'
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 

plot_transMat_heatmap_from_adata(split1, split2, total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)

## utv
method = 'utv'
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = scv.read(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 

plot_transMat_heatmap_from_adata(split1, split2, total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)

## sct
method = 'sct'
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 

plot_transMat_heatmap_from_adata(split1, split2, total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)

## velovi
method = 'velovi'
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split1_outputAdded_v2.h5ad")
split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split2_outputAdded_v2.h5ad")
total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_total_outputAdded_v2.h5ad")

plot_transMat_heatmap_from_adata(split1, split2, total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)

## velovi woprep
method = 'velovi'
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/wopreprocess/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_outputAdded_v2.h5ad")
split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_outputAdded_v2.h5ad")
total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_outputAdded_v2.h5ad")

plot_transMat_heatmap_from_adata(split1, split2, total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)
