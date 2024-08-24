# unitvelo virtual environment

import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import datetime
import random
import numpy as np

sct_seed = 615

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1
velo_config.GPU = -1

os.environ["TF_USE_LEGACY_KERAS"]="1"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

adata = sc.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup9_larry-subset-2block_split1.h5ad")

random.seed(sct_seed)
np.random.seed(sct_seed)

adata = utv.run_model(adata,
                      'state_info', 
                      config_file=velo_config)

print_message_with_time("########### Saving")
adata.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup9_unitvelo_larry-subset-2block_split1.h5ad")
print_message_with_time("########### Total data wrote")

