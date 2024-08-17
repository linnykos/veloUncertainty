import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sctour as sct

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from sctour_misc import *

method = 'sct'
dataset_long = 'pancreas'
dataset_short = 'pan'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_seed317_split1_v2.h5ad') # 

tnode_split1 = sct.predict.load_model(data_folder+'v2_pancreas/sct/tnode_pan_sct_seed317_split1_v2.pth')
velo_mat_split1 = compute_sctour_velocity(tnode_split1, timestep=1/100)

np.sum(split1.layers['velocity']+velo_mat_split1) # 0.0
split1.layers['velocity'] = velo_mat_split1
split1.write(data_folder+'v2_pancreas/sct/adata_pan_sct_seed317_split1_v2.h5ad')
split1.write(data_folder+'v2_pancreas/sct/adata_pan_sct_split1_v2.h5ad')

