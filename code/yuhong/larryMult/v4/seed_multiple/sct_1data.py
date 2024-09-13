import sctour as sct
import scanpy as sc
import numpy as np
import torch
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *


def sct_larryMult(split_seed):
    method = 'sct'
    dataset_long = 'larryMult'
    dataset_short = 'larryMult'
    sct_seed = 615
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    split_version = 'split1'
    read_data_and_run_sct(dataset_short,dataset_long,method,data_folder,savedata_folder,split_version,split_seed,sct_seed)
    split_version = 'split2'
    read_data_and_run_sct(dataset_short,dataset_long,method,data_folder,savedata_folder,split_version,split_seed,sct_seed)
    print('All done with seed='+str(split_seed))

sct_larryMult(320)
sct_larryMult(323)
sct_larryMult(326)
sct_larryMult(329)
