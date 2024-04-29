import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

label='celltype'
adata = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utv_preprocess.h5ad')

total_res = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utvgenes_total.h5ad")
s1_res317 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utvgenes_seed317_split1.h5ad")
s2_res317 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utvgenes_seed317_split2.h5ad")

common_genes_317 = np.intersect1d(s1_res317.var['features'], s2_res317.var['features'])
Ngenes_317s1 = len(s1_res317.var['features']) # 1425
Ngenes_317s2 = len(s2_res317.var['features']) # 1442
Ngenes_317common = len(common_genes_317) # 1331

df1_317 = pd.DataFrame(s1_res317.layers['velocity'], columns=s1_res317.var['features'].tolist())
df2_317 = pd.DataFrame(s2_res317.layers['velocity'], columns=s2_res317.var['features'].tolist())
cos_sim_seed317 = np.diag(cosine_similarity(df1_317[common_genes_317],df2_317[common_genes_317]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed317 = np.mean(cos_sim_seed317)
plt.axvline(mean_seed317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-.75, 3000, 'mean cosine similarity = '+str(mean_seed317), color='blue', fontsize=10)
plt.text(-.75, 2000, 'split1 number of genes = '+str(Ngenes_317s1), color='blue', fontsize=10)
plt.text(-.75, 1000, 'split2 number of genes = '+str(Ngenes_317s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, ery+utv_utvgenes, Ngenes='+str(Ngenes_317common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed317_cos_similarity_hist.png')
plt.clf()

total_res.obs['cos_sim_317'] = cos_sim_seed317
total_res.obs['cos_sim_317'] = pd.DataFrame(total_res.obs['cos_sim_317'])

total_res.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(total_res,color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed317_cos_similarity_preumap.png")
del total_res.obsm['X_umap']
sc.tl.umap(total_res)
scv.pl.velocity_embedding_stream(total_res,color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed317_cos_similarity_uncorrected.png")
scv.pl.velocity_embedding_stream(total_res,color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed317total_uncorrected.png")
del total_res.obsm['X_umap']



###########################################################
s1_res320 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utvgenes_seed320_split1.h5ad")
s2_res320 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utvgenes_seed320_split2.h5ad")

common_genes_320 = np.intersect1d(s1_res320.var['features'], s2_res320.var['features'])
Ngenes_320s1 = len(s1_res320.var['features']) # 1433
Ngenes_320s2 = len(s2_res320.var['features']) # 1446
Ngenes_320common = len(common_genes_320) # 1342

df1_320 = pd.DataFrame(s1_res320.layers['velocity'], columns=s1_res320.var['features'].tolist())
df2_320 = pd.DataFrame(s2_res320.layers['velocity'], columns=s2_res320.var['features'].tolist())
cos_sim_seed320 = np.diag(cosine_similarity(df1_320[common_genes_320],df2_320[common_genes_320]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed320, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed320 = np.mean(cos_sim_seed320)
plt.axvline(mean_seed320, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-.75, 3000, 'mean cosine similarity = '+str(mean_seed320), color='blue', fontsize=10)
plt.text(-.75, 2000, 'split1 number of genes = '+str(Ngenes_320s1), color='blue', fontsize=10)
plt.text(-.75, 1000, 'split2 number of genes = '+str(Ngenes_320s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, ery+utv_utvgenes, Ngenes='+str(Ngenes_320common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed320_cos_similarity_hist.png')
plt.clf()

total_res.obs['cos_sim_320'] = cos_sim_seed320
total_res.obs['cos_sim_320'] = pd.DataFrame(total_res.obs['cos_sim_320'])

total_res.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(total_res,color="cos_sim_320",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed320_cos_similarity_preumap.png")
del total_res.obsm['X_umap']
sc.tl.umap(total_res)
scv.pl.velocity_embedding_stream(total_res,color="cos_sim_320",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed320_cos_similarity_uncorrected.png")
scv.pl.velocity_embedding_stream(total_res,color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_seed320total_uncorrected.png")
del total_res.obsm['X_umap']