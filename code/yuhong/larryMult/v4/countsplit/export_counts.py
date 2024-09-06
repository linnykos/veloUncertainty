import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import scvelo as scv

import scipy.sparse as sp
from scipy.io import mmwrite

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import read_raw_adata,print_message_with_time

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"


## larryMult
dataset_short = 'larryMult'
dataset_long = 'larryMult'

total = read_raw_adata(dataset_short) # 2315 × 23420
total_S = total.layers['spliced'].copy().astype(float)
total_U = total.layers['unspliced'].copy().astype(float)

# Save the sparse matrix to a Matrix Market file
mmwrite(data_folder+'v4_larryMult/larryMult_spliced.mtx', total_S)
mmwrite(data_folder+'v4_larryMult/larryMult_unspliced.mtx', total_U)

#############################
## erythroid
total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
total_S = total.layers['spliced'].copy().astype(float)
total_U = total.layers['unspliced'].copy().astype(float)

mmwrite(data_folder+'v2_erythroid/ery_spliced.mtx', total_S)
mmwrite(data_folder+'v2_erythroid/ery_unspliced.mtx', total_U)

#############################
## pancreas
total = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
total_S = total.layers['spliced'].copy().astype(float)
total_U = total.layers['unspliced'].copy().astype(float)

mmwrite(data_folder+'v2_pancreas/pan_spliced.mtx', total_S)
mmwrite(data_folder+'v2_pancreas/pan_unspliced.mtx', total_U)


