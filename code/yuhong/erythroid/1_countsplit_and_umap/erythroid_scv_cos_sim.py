import scvelo as scv
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bbknn

adata_split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_split1_seurat.h5ad')
adata_split2_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_split2_seurat.h5ad')
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed317_total_seurat.h5ad')

scv.pp.normalize_per_cell(adata_split1_seed317)
scv.pp.log1p(adata_split1_seed317)
scv.pp.moments(adata_split1_seed317, n_pcs=30, n_neighbors=30)
## batch correction
bbknn.bbknn(adata_split1_seed317, batch_key='sequencing.batch')
adata_split1_seed317.X = adata_split1_seed317.X.toarray()
bbknn.ridge_regression(adata_split1_seed317, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split1_seed317)
bbknn.bbknn(adata_split1_seed317, batch_key='sequencing.batch')
print("Batch correction for seed317 split1 done!")
sc.pp.neighbors(adata_split1_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_seed317)
scv.tl.recover_dynamics(adata_split1_seed317)
scv.tl.velocity(adata_split1_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split1_seed317)
scv.pl.velocity_embedding_stream(adata_split1_seed317, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed317_split1.png")
print("**************** seed317 split1 processed! ****************")


scv.pp.normalize_per_cell(adata_split2_seed317)
scv.pp.log1p(adata_split2_seed317)
scv.pp.moments(adata_split2_seed317, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_split2_seed317, batch_key='sequencing.batch')
adata_split2_seed317.X = adata_split2_seed317.X.toarray()
bbknn.ridge_regression(adata_split2_seed317, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split2_seed317)
bbknn.bbknn(adata_split2_seed317, batch_key='sequencing.batch')
print("Batch correction done for seed317 split2 total counts!")
sc.pp.neighbors(adata_split2_seed317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_seed317)
scv.tl.recover_dynamics(adata_split2_seed317)
scv.tl.velocity(adata_split2_seed317, mode="dynamical")
scv.tl.velocity_graph(adata_split2_seed317)
scv.pl.velocity_embedding_stream(adata_split2_seed317, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed317_split2.png")
print("**************** seed317 split2 processed! ****************")

# replace nan's to 0's in layers['velocity']
adata_split1_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed317.layers['velocity'], nan=0)
adata_split2_seed317.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed317.layers['velocity'], nan=0)
Ngenes_317s1 = np.sum(~np.isnan(adata_split1_seed317.layers['velocity'][0]))
Ngenes_317s2 = np.sum(~np.isnan(adata_split2_seed317.layers['velocity'][0]))
Ngenes_317common = np.sum(np.isnan(adata_split1_seed317.layers["velocity"][0] + adata_split2_seed317.layers["velocity"][0])==0)
# cosine similarity
cos_sim_seed317 = np.diag(cosine_similarity(adata_split1_seed317.layers["velocity_rmNA"],
                                            adata_split2_seed317.layers["velocity_rmNA"]))
print("**************** seed317 erythroid cosine similarity computed! ****************")

# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed317 = np.mean(cos_sim_seed317)
plt.axvline(mean_seed317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-1, 320, 'mean cosine similarity = '+str(mean_seed317), color='blue', fontsize=10)
plt.text(-1, 300, 'split1 number of genes = '+str(Ngenes_317s1), color='blue', fontsize=10)
plt.text(-1, 280, 'split2 number of genes = '+str(Ngenes_317s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, ery+scv, Ngenes='+str(Ngenes_317common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed317_cos_similarity_hist.png')
plt.clf()

# total counts data process
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
scv.tl.recover_dynamics(adata_total)
scv.tl.velocity(adata_total, mode="dynamical")
scv.tl.velocity_graph(adata_total)
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed317_total.png")
print("**************** total counts processed! ****************")

# add cosine similarities to total counts object
adata_total.obs["cos_sim_seed317"] = cos_sim_seed317
adata_total.obs["cos_sim_seed317"] = pd.DataFrame(adata_total.obs["cos_sim_seed317"])
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="cos_sim_seed317",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed317_cos_similarity.png")
print("**************** seed317 cosine similarity plotted! ****************")

# read split counts data
adata_split1_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed320_split1_seurat.h5ad')
adata_split2_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed320_split2_seurat.h5ad')
print("**************** read seed320 split counts! ****************")

## process split1
scv.pp.normalize_per_cell(adata_split1_seed320)
scv.pp.log1p(adata_split1_seed320)
scv.pp.moments(adata_split1_seed320, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_split1_seed320, batch_key='sequencing.batch')
adata_split1_seed320.X = adata_split1_seed320.X.toarray()
bbknn.ridge_regression(adata_split1_seed320, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split1_seed320)
bbknn.bbknn(adata_split1_seed320, batch_key='sequencing.batch')
print("Batch correction done for seed320 split1 counts!")
sc.pp.neighbors(adata_split1_seed320, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_seed320)
scv.tl.recover_dynamics(adata_split1_seed320)
scv.tl.velocity(adata_split1_seed320, mode="dynamical")
scv.tl.velocity_graph(adata_split1_seed320)
scv.pl.velocity_embedding_stream(adata_split1_seed320, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed320_split1_2000.png")
print("**************** seed320 split1 processed! ****************")

## process split2
scv.pp.normalize_per_cell(adata_split2_seed320)
scv.pp.log1p(adata_split2_seed320)
scv.pp.moments(adata_split2_seed320, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_split2_seed320, batch_key='sequencing.batch')
adata_split2_seed320.X = adata_split2_seed320.X.toarray()
bbknn.ridge_regression(adata_split2_seed320, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split2_seed320)
bbknn.bbknn(adata_split2_seed320, batch_key='sequencing.batch')
print("Batch correction done for seed320 split2 counts!")
sc.pp.neighbors(adata_split2_seed320, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_seed320)
scv.tl.recover_dynamics(adata_split2_seed320)
scv.tl.velocity(adata_split2_seed320, mode="dynamical")
scv.tl.velocity_graph(adata_split2_seed320)
scv.pl.velocity_embedding_stream(adata_split2_seed320, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed320_split2_2000.png")
print("**************** seed320 split2 processed! ****************")

# replace nan's to 0's in layers['velocity']
adata_split1_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split1_seed320.layers['velocity'], nan=0)
adata_split2_seed320.layers["velocity_rmNA"] = np.nan_to_num(adata_split2_seed320.layers['velocity'], nan=0)

Ngenes_320s1 = np.sum(~np.isnan(adata_split1_seed320.layers['velocity'][0]))
Ngenes_320s2 = np.sum(~np.isnan(adata_split2_seed320.layers['velocity'][0]))
Ngenes_320common = np.sum(np.isnan(adata_split1_seed320.layers["velocity"][0] + adata_split2_seed320.layers["velocity"][0])==0)
# cosine similarity
cos_sim_seed320 = np.diag(cosine_similarity(adata_split1_seed320.layers["velocity_rmNA"],
                                            adata_split2_seed320.layers["velocity_rmNA"]))
print("**************** seed320 erythroid cosine similarity computed! ****************")

# Create histogram
plt.clf()
plt.hist(cos_sim_seed320, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed320 = np.mean(cos_sim_seed320)
plt.axvline(mean_seed320, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-1, 320, 'mean cosine similarity = '+str(mean_seed320), color='blue', fontsize=10)
plt.text(-1, 300, 'split1 number of genes = '+str(Ngenes_320s1), color='blue', fontsize=10)
plt.text(-1, 280, 'split2 number of genes = '+str(Ngenes_320s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, ery+scv, Ngenes='+str(Ngenes_320common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed320_cos_similarity_hist.png')
plt.clf()

# add cosine similarities to total counts object
adata_total.obs["cos_sim_seed320"] = cos_sim_seed320
adata_total.obs["cos_sim_seed320"] = pd.DataFrame(adata_total.obs["cos_sim_seed320"])
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="cos_sim_seed320",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/scvelo_seed320_cos_similarity_2000.png")
print("**************** seed320 cosine similarity plotted! ****************")


