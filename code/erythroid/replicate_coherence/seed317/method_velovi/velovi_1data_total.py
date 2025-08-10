data_version = "total"
dataset_short = 'ery'
dataset_long = 'erythroid'
method = 'velovi'
split_seed = 317
celltype_label = 'celltype'
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

velovi_run_model(data_version=data_version,dataset_short=dataset_short,dataset_long=dataset_long,method=method,data_folder=data_folder,split_seed=split_seed)
