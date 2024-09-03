split_seed = 317
split_version = 'split2'
method = 'sct'
dataset_long = 'erythroid'
dataset_short = 'ery'
sct_seed = 615

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'#+'seed'+str(split_seed)+'/'
savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'

import sctour as sct
import scanpy as sc
import numpy as np
import torch
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *

######################################################

read_data_and_run_sct(dataset_short,dataset_long,method,data_folder,savedata_folder,split_version,split_seed,sct_seed)