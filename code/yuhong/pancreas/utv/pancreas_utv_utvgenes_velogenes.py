import scvelo as scv
import numpy as np
import pandas as pd

## seed=317
adata_total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_total_seurat.h5ad')
#adata_total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utv_preprocess.h5ad')
#adata_total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_total.h5ad') # 3696 × 1945
split1_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_split1.h5ad') # 3696 × 734
split2_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_split2.h5ad')

common_genes_317 = np.intersect1d(split1_seed317_res.var['features'], split2_seed317_res.var['features'])
indices_317 = [i for i, value in enumerate(adata_total_res.var['features']) if value in common_genes_317]

np.savetxt('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup9_corr/pan_utv_genes_seed317.csv', 
           indices_317, delimiter=',')

## seed=320
split1_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split1.h5ad') # 3696 × 734
split2_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split2.h5ad')

common_genes_320 = np.intersect1d(split1_seed320_res.var['features'], split2_seed320_res.var['features'])
indices_320 = [i for i, value in enumerate(adata_total_res.var['features']) if value in common_genes_320]

np.savetxt('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup9_corr/pan_utv_genes_seed320.csv', 
           indices_320, delimiter=',')


