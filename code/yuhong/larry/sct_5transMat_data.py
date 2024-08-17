import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity
import sctour as sct

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from sctour_misc import *

method = 'sct'
dataset_long = 'larry'
dataset_short = 'larry'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

tnode_total = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_total_v2.pth')
tnode_split1 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split1_v2.pth')
tnode_split2 = sct.predict.load_model(data_folder+'v2_'+dataset_long+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split2_v2.pth')

total = tnode_total.adata
split1 = tnode_split1.adata
split2 = tnode_split2.adata

total.layers['velocity'] = compute_sctour_velocity(tnode_total, timestep=0.42)
split1.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=0.42)
split2.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=0.42)

total.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_split2_v2.h5ad')
