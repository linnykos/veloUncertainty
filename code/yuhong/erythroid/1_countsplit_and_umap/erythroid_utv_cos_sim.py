import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad")
#adata = scv.datasets.gastrulation_erythroid()
scv.pp.filter_genes(adata, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) 
### Extracted 2000 highly variable genes.
scv.pp.log1p(adata)

total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_total.h5ad')
s1_res317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed317_split1.h5ad')
s2_res317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed317_split2.h5ad')

common_genes_317 = np.intersect1d(s1_res317.var['features'], s2_res317.var['features']) # 1331
# indices_seed317 = s1_res317.var['features'][s1_res317.var['features'].isin(set(common_genes_317))].index.tolist()
Ngenes_317s1 = len(s1_res317.var['features'])
Ngenes_317s2 = len(s2_res317.var['features'])
Ngenes_317common = len(common_genes_317)

df1_317 = pd.DataFrame(s1_res317.layers['velocity'], columns=s1_res317.var['features'].tolist())
df2_317 = pd.DataFrame(s2_res317.layers['velocity'], columns=s2_res317.var['features'].tolist())
cos_sim_seed317 = np.diag(cosine_similarity(df1_317[common_genes_317],df2_317[common_genes_317]))
# cos_sim_seed317 = np.diag(cosine_similarity(df1_317[indices_seed317],df2_317[indices_seed317]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed317 = np.mean(cos_sim_seed317)
plt.axvline(mean_seed317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-.75, 525, 'mean cosine similarity = '+str(mean_seed317), color='blue', fontsize=10)
plt.text(-.75, 500, 'split1 number of genes = '+str(Ngenes_317s1), color='blue', fontsize=10)
plt.text(-.75, 475, 'split2 number of genes = '+str(Ngenes_317s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, ery+utv, Ngenes='+str(Ngenes_317common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed317_cos_similarity_hist.png')
plt.clf()

total_res.obs['cos_sim_317'] = cos_sim_seed317
total_res.obs['cos_sim_317'] = pd.DataFrame(total_res.obs['cos_sim_317'])
scv.pl.velocity_embedding_stream(total_res,color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed317_cos_similarity_nonumap.png")
scv.pl.velocity_embedding_stream(total_res,color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_uncorrected_nonumap.png")

total_res.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(total_res,basi='umap',color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed317_cos_similarity_preumap.png")
print("seed317 cosine similarity done!")


#####################################################################
s1_res320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed320_split1.h5ad')
s2_res320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed320_split2.h5ad')
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_total_seurat.h5ad')

common_genes_320 = np.intersect1d(s1_res320.var['features'], s2_res320.var['features']) # 1342
# indices_seed317 = s1_res317.var['features'][s1_res317.var['features'].isin(set(common_genes_317))].index.tolist()
Ngenes_320s1 = len(s1_res320.var['features']) # 1433
Ngenes_320s2 = len(s2_res320.var['features']) # 1446
Ngenes_320common = len(common_genes_320)

df1_320 = pd.DataFrame(s1_res320.layers['velocity'], columns=s1_res320.var['features'].tolist())
df2_320 = pd.DataFrame(s2_res320.layers['velocity'], columns=s2_res320.var['features'].tolist())
cos_sim_seed320 = np.diag(cosine_similarity(df1_320[common_genes_320],df2_320[common_genes_320]))
# cos_sim_seed317 = np.diag(cosine_similarity(df1_317[indices_seed317],df2_317[indices_seed317]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed320, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed320 = np.mean(cos_sim_seed320)
plt.axvline(mean_seed320, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-.75, 525, 'mean cosine similarity = '+str(mean_seed320), color='blue', fontsize=10)
plt.text(-.75, 500, 'split1 number of genes = '+str(Ngenes_320s1), color='blue', fontsize=10)
plt.text(-.75, 475, 'split2 number of genes = '+str(Ngenes_320s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, ery+utv, Ngenes='+str(Ngenes_320common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed320_cos_similarity_hist.png')
plt.clf()

total_res.obs['cos_sim_320'] = cos_sim_seed320
total_res.obs['cos_sim_320'] = pd.DataFrame(total_res.obs['cos_sim_320'])
scv.pl.velocity_embedding_stream(total_res,color="cos_sim_320",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed320_cos_similarity_nonumap.png")

total_res.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(total_res,basi='umap',color="cos_sim_320",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed320_cos_similarity_preumap.png")

print("seed320 cosine similarity done!")


## plot cosine similarity on umap from adata processed by scv and corrected by bbknn
import bbknn
import scanpy as sc

adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_total_seurat.h5ad')

scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_total, batch_key='sequencing.batch')
adata_total.X = adata_total.X.toarray()
bbknn.ridge_regression(adata_total, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_total)
bbknn.bbknn(adata_total, batch_key='sequencing.batch')
print("Batch correction done for total counts!")
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
adata_total.obs['cos_sim_317'] = cos_sim_seed317
adata_total.obs['cos_sim_320'] = cos_sim_seed320

total_res_317 = total_res.copy()
total_res_317.obsm['X_umap'] = adata_total.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(total_res_317, basis='umap',color="cos_sim_317",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed317_cos_similarity_scvumapbbknn.png")
total_res_320 = total_res.copy()
total_res_320.obsm['X_umap'] = adata_total.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(total_res_320, basis='umap',color="cos_sim_320",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_seed320_cos_similarity_scvumapbbknn.png")

