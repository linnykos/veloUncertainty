import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_velovi import *

split_seed = 317
method = 'velovi_woprep'
dataset_short = 'pan'
dataset_long = 'pancreas'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

## intrinsic uncertainty
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
vae_total = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_total_v4.pt', total)

print_message_with_time("############## Plot intrinsic uncertainty for total")
compute_intrinisic_uncertainty(adata_in=total,vae=vae_total,dataset=dataset_short,fig_folder=fig_folder,data_version='total',sample_seed=2216,n_samples=100,recompute=True,celltype_label=None)

#######################################
## extrinsic uncertainty
print_message_with_time("############## Plot extrinsic uncertainty for total")
compute_extrinisic_uncertainty(adata_in=total,vae=vae_total,dataset=dataset_short,fig_folder=fig_folder,data_version='total',sample_seed=2216,n_samples=25,recompute=True,celltype_label=None)

total.write(data_folder+'v4_pancreas/seed317/velovi_woprep/adata_pan_velovi_woprep_total_v4_outputAdded_uncertaintyAdded.h5ad')

