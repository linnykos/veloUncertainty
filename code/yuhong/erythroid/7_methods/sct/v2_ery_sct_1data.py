import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from sctour_misc import *

sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/sct/" 

print_message_with_time("########### Start to read total, split1 and split2")
total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
adata_split1 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split1_allgenes.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split2_allgenes.h5ad')
adata_split1.var['Accession'] = total.var['Accession'].copy()
adata_split1.var['Accession'] = total.var['Accession'].copy()

def train_sct_and_return_tnode(adata, sct_seed=615):
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    adata.X = adata.X.astype(np.float32)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)
    tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
    return tnode

print_message_with_time("########### Start to train model for split1")
tnode_split1 = train_sct_and_return_tnode(adata_split1)
print_message_with_time("########### Start to compute velocity for split1")
diff_mat_split1 = compute_sctour_velocity(tnode_split1, timestep=1/100)
print_message_with_time("########### Split1 velocity computed, start to write data")
adata_split1.layers['velocity'] = diff_mat_split1
adata_split1.write(data_folder+'v2_erythroid/sct/adata_ery_sct_seed317_split1_v2.h5ad')
tnode_split1.save_model(save_dir=data_folder+'v2_erythroid/sct/', save_prefix='tnode_ery_sct_seed317_split1_v2')
print_message_with_time("########### Split1 data wrote")

print_message_with_time("########### Start to train model for split2")
tnode_split2 = train_sct_and_return_tnode(adata_split2)
print_message_with_time("########### Start to compute velocity for split2")
diff_mat_split2 = compute_sctour_velocity(tnode_split2, timestep=1/100)
print_message_with_time("########### Split2 velocity computed, start to write data")
adata_split2.layers['velocity'] = diff_mat_split2
adata_split2.write(data_folder+'v2_erythroid/sct/adata_ery_sct_seed317_split2_v2.h5ad')
tnode_split2.save_model(save_dir=data_folder+'v2_erythroid/sct/', save_prefix='tnode_ery_sct_seed317_split2_v2')
print_message_with_time("########### Split2 data wrote")

print_message_with_time("########### Start to train model for total")
tnode_total = train_sct_and_return_tnode(total)
print_message_with_time("########### Start to compute velocity for split1")
diff_mat_total = compute_sctour_velocity(tnode_total, timestep=1/100)
print_message_with_time("########### Total velocity computed, start to write data")
total.layers['velocity'] = diff_mat_total
total.write(data_folder+'v2_erythroid/sct/adata_ery_sct_total_v2.h5ad')
tnode_total.save_model(save_dir=data_folder+'v2_erythroid/sct/', save_prefix='tnode_ery_sct_total_v2')
print_message_with_time("########### Total data wrote")

