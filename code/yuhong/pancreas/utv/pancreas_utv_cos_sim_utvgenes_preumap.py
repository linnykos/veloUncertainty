import unitvelo as utv
import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata=scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utv_preprocess.h5ad')

## seed317
split1_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_split1.h5ad') # 3696 × 734
split2_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_split2.h5ad') # 3696 × 731

common_genes_317 = np.intersect1d(split1_seed317_res.var['features'], split2_seed317_res.var['features']) # 649
Ngenes_317s1 = len(split1_seed317_res.var['features'])
Ngenes_317s2 = len(split2_seed317_res.var['features'])
Ngenes_317common = len(common_genes_317)

df1_317 = pd.DataFrame(split1_seed317_res.layers['velocity'], columns=split1_seed317_res.var['features'].tolist())
df2_317 = pd.DataFrame(split2_seed317_res.layers['velocity'], columns=split2_seed317_res.var['features'].tolist())
cos_sim_seed317 = np.diag(cosine_similarity(df1_317[common_genes_317],df2_317[common_genes_317]))

adata.obs['cos_sim_317'] = cos_sim_seed317
adata.obs['cos_sim_317'] = pd.DataFrame(adata.obs['cos_sim_317'])
scv.pl.velocity_embedding_stream(adata,basis="umap",color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/pan_utv_seed317_cos_sim_preumap_utvgenes.png")

## seed320
split1_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split1.h5ad') 
split2_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split2.h5ad') 

common_genes_320 = np.intersect1d(split1_seed320_res.var['features'], split2_seed320_res.var['features']) # 649
Ngenes_320s1 = len(split1_seed320_res.var['features'])
Ngenes_320s2 = len(split2_seed320_res.var['features'])
Ngenes_320common = len(common_genes_320)

df1_320 = pd.DataFrame(split1_seed320_res.layers['velocity'], columns=split1_seed320_res.var['features'].tolist())
df2_320 = pd.DataFrame(split2_seed320_res.layers['velocity'], columns=split2_seed320_res.var['features'].tolist())
cos_sim_seed320 = np.diag(cosine_similarity(df1_320[common_genes_320],df2_320[common_genes_320]))

adata.obs['cos_sim_320'] = cos_sim_seed320
adata.obs['cos_sim_320'] = pd.DataFrame(adata.obs['cos_sim_320'])
scv.pl.velocity_embedding_stream(adata,basis="umap",color="cos_sim_320",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/pan_utv_seed320_cos_sim_preumap_utvgenes.png")

