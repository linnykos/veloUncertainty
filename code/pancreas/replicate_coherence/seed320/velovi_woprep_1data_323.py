split_seed=323
method = 'velovi_woprep'
dataset_short = 'pan'
dataset_long = 'pancreas'
celltype_label = 'clusters'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_velovi import *

data_version = 'split1'
velovi_run_model(data_version=data_version,dataset_long=dataset_long,dataset_short=dataset_short,method=method,data_folder=data_folder,split_seed=split_seed)
data_version = 'split2'
velovi_run_model(data_version=data_version,dataset_long=dataset_long,dataset_short=dataset_short,method=method,data_folder=data_folder,split_seed=split_seed)
