import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity

def train_sct(adata):
    adata.X = adata.X.astype(np.float32)
    adata.layers['spliced'] = adata.layers['spliced'].astype(np.float32)
    adata.layers['unspliced'] = adata.layers['unspliced'].astype(np.float32)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)
    tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    return adata

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_sct/"
total = sc.read(data_folder+"ery_sct_preprocess.h5ad")

s1_317 = sc.read(data_folder+"erythroid_seed317_split1_seurat.h5ad")
s1_317_res = train_sct(s1_317)
s1_317_res.write(filename=data_folder+"erythroid_seed317_split1.h5ad")

s2_317 = sc.read(data_folder+"erythroid_seed317_split2_seurat.h5ad")
s2_317_res = train_sct(s2_317)
s2_317_res.write(filename=data_folder+"erythroid_seed317_split2.h5ad")

s1_320 = sc.read(data_folder+"erythroid_seed320_split1_seurat.h5ad")
s1_320_res = train_sct(s1_320)
s1_320_res.write(filename=data_folder+"erythroid_seed320_split1.h5ad")

s2_320 = sc.read(data_folder+"erythroid_seed320_split2_seurat.h5ad")
s2_320_res = train_sct(s2_320)
s2_320_res.write(filename=data_folder+"erythroid_seed320_split2.h5ad")



