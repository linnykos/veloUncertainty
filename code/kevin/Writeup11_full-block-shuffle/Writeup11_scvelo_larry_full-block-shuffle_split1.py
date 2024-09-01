import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
import datetime

sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

adata = sc.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/Writeup11_larry_full-block_split1.h5ad")

torch.manual_seed(sct_seed)
random.seed(sct_seed)
np.random.seed(sct_seed)

scv.pp.filter_and_normalize(adata, 
                            min_shared_counts=20, 
                            n_top_genes=2000)
scv.pp.moments(adata, 
               n_pcs=30, 
               n_neighbors=30)
sc.tl.pca(adata, 
          svd_solver="arpack")
sc.pp.neighbors(adata, 
                n_neighbors=15, 
                n_pcs=40) # used to be n_neighbors=10
sc.tl.umap(adata)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, 
                mode="dynamical")
scv.tl.velocity_graph(adata)

adata.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/Writeup11_scvelo_larry_full-block_split1.h5ad")
print_message_with_time("########### Split1 data wrote")

