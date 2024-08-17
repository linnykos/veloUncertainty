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

method = 'scv'
dataset_long = 'erythroid'
dataset_short = 'ery'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_total_v2.h5ad')
split1 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_split1_v2.h5ad')
split2 = scv.read(data_folder+'v2_erythroid/scv/adata_ery_scv_split2_v2.h5ad')

plot_transMat_heatmap_from_adata(split1, split2, total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)
