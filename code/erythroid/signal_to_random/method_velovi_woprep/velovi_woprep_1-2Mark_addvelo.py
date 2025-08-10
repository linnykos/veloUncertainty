method_prefix = 'velovi_woprep'
gene_set_name = 'Mark'

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_velovi import *
from v4_functions import *

def addvelo_velovi_ery_Mark(method_prefix, gene_set_name, split_seed, total_velo=False):
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method = method_prefix + '_' + gene_set_name
    data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    ## read data
    split1 = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split1'+'.h5ad')
    split2 = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split2'+'.h5ad')
    vae_split1 = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'split1.pt', split1)
    vae_split2 = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'split2.pt', split2)
    if total_velo: 
        total = sc.read_h5ad(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'total'+'.h5ad')
        vae_total = VELOVI.load(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+'total.pt', total)
    ## add velovi outputs to adata
    print_message_with_time("############## Add velovi outputs to adata")
    add_velovi_outputs_to_adata(split1, vae_split1)
    add_velovi_outputs_to_adata(split2, vae_split2)
    if total_velo: add_velovi_outputs_to_adata(total, vae_total)
    ## compute umap
    print_message_with_time("############## Compute umap")
    if total_velo: compute_umap_ery(total)
    compute_umap_ery(split1)
    compute_umap_ery(split2)
    split1.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split1'+'_outputAdded.h5ad')
    split2.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'split2'+'_outputAdded.h5ad')
    if total_velo: total.write(savedata_folder+'adata_'+dataset_short+'_'+method+'_'+'total'+'_outputAdded.h5ad')
    print('####################### seed'+str(split_seed)+' done.')

total_velo = True
for split_seed in [317,320,323,326,329]:
    addvelo_velovi_ery_Mark(method_prefix=method_prefix, gene_set_name=gene_set_name, split_seed=split_seed, total_velo=total_velo)
    total_velo = False

print('####################### All done.')
