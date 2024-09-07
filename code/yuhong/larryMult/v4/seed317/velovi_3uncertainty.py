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
method = 'velovi'
dataset_short = 'larryMult'
dataset_long = 'larryMult'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)

vae_split1 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split1_v4.pt', split1)
vae_split2 = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_split2_v4.pt', split2)
vae_total = VELOVI.load(data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/vae_'+dataset_short+'_'+method+'_total_v4.pt', total)

######################################################
## intrinsic uncertainty
compute_intrinisic_uncertainty(adata_in=split1,vae=vae_split1,dataset=dataset_short,fig_folder=fig_folder,data_version='split1')
compute_intrinisic_uncertainty(adata_in=split2,vae=vae_split2,dataset=dataset_short,fig_folder=fig_folder,data_version='split2')
compute_intrinisic_uncertainty(adata_in=total,vae=vae_total,dataset=dataset_short,fig_folder=fig_folder,data_version='total')

## extrinsic uncertainty
compute_extrinisic_uncertainty(adata_in=split1,vae=vae_split1,dataset=dataset_short,fig_folder=fig_folder,data_version='split1')
compute_extrinisic_uncertainty(adata_in=split2,vae=vae_split2,dataset=dataset_short,fig_folder=fig_folder,data_version='split2')
compute_extrinisic_uncertainty(adata_in=total,vae=vae_total,dataset=dataset_short,fig_folder=fig_folder,data_version='total')

## permutation score
compute_permutation_score(adata=split1,vae=vae_split1,dataset=dataset_short,fig_folder=fig_folder,data_version='split1')
compute_permutation_score(adata=split2,vae=vae_split2,dataset=dataset_short,fig_folder=fig_folder,data_version='split2')
compute_permutation_score(adata=total,vae=vae_total,dataset=dataset_short,fig_folder=fig_folder,data_version='total')