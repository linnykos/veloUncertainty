split_seed = 317
dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'utv'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
data_version = 'split2'
celltype_label = 'cluster_name'

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_utv import *

# the below script uses the environment: "utvClone"
# larry/erythroid configuration
velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True # default value
# (bool) linear regression $R^2$ on extreme quantile (default) or full data (adjusted)
# valid when self.VGENES = 'basic'
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1  # default initialization
# (float) threshold of R2 at later stage of the optimization proces
# to capture the dynamics of more genes beside initially selected velocity genes
# self.AGENES_R2 = 1 will switch to origianl mode with no amplification
os.environ["TF_USE_LEGACY_KERAS"]="1"

utv_run_model_v4(data_version,dataset_short,method,velo_config,data_folder,split_seed,celltype_label)
