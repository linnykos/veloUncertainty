import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_transMat import *

method = 'velovi'
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/wopreprocess/"

split1=sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_outputAdded_v2.h5ad")
split2=sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_outputAdded_v2.h5ad")
total=sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_outputAdded_v2.h5ad")

plot_transMat_heatmap_from_adata(split1, split2, total, method=method, dataset_short=dataset_short, fig_folder=fig_folder)
