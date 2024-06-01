##### cosine similarity
import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata_total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes2000_total.h5ad') # 3696 × 1945
split1_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes2000_seed317_split1.h5ad') # 3696 × 734
split2_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes2000_seed317_split2.h5ad') # 3696 × 731

common_genes_317 = np.intersect1d(split1_seed317_res.var['features'], split2_seed317_res.var['features']) # 649
Ngenes_317s1 = len(split1_seed317_res.var['features'])
Ngenes_317s2 = len(split2_seed317_res.var['features'])
Ngenes_317common = len(common_genes_317)

# Here assuming order of genes in layers['velocity'] is the same as in var['features']
df1_317 = pd.DataFrame(split1_seed317_res.layers['velocity'], columns=split1_seed317_res.var['features'].tolist())
df2_317 = pd.DataFrame(split2_seed317_res.layers['velocity'], columns=split2_seed317_res.var['features'].tolist())
cos_sim_seed317 = np.diag(cosine_similarity(df1_317[common_genes_317],df2_317[common_genes_317]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed317 = np.mean(cos_sim_seed317)
plt.axvline(mean_seed317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.25, 500, 'mean cosine similarity = '+str(mean_seed317), color='blue', fontsize=10)
plt.text(.25, 450, 'split1 number of genes = '+str(Ngenes_317s1), color='blue', fontsize=10)
plt.text(.25, 400, 'split2 number of genes = '+str(Ngenes_317s2), color='blue', fontsize=10)
plt.text(.25, 350, 'velocities computed with VGENES specified (highly_variable=T)', color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+utv_utvgenes2000, Ngenes='+str(Ngenes_317common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/utvgenes2000_seed317_cos_sim_hist.png')
plt.clf()

## Plot on UMAP
adata_total_res.obs['cos_sim_317'] = cos_sim_seed317
adata_total_res.obs['cos_sim_317'] = pd.DataFrame(adata_total_res.obs['cos_sim_317'])
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/utvgenes2000_seed317_cos_sim_preumap.png")

print("seed317 cosine similarity done!")
#############################
adata_total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes2000_total.h5ad') # 3696 × 1945
split1_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes2000_seed320_split1.h5ad')
split2_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes2000_seed320_split2.h5ad')

Ngenes_320s1 = len(split1_seed320_res.var['features']) # 741
Ngenes_320s2 = len(split2_seed320_res.var['features']) # 744
common_genes_320 = np.intersect1d(split1_seed320_res.var['features'], split2_seed320_res.var['features'])
Ngenes_320common = len(common_genes_320) # 654

df1_320 = pd.DataFrame(split1_seed320_res.layers['velocity'], columns=split1_seed320_res.var['features'].tolist())
df2_320 = pd.DataFrame(split2_seed320_res.layers['velocity'], columns=split2_seed320_res.var['features'].tolist())
cos_sim_seed320 = np.diag(cosine_similarity(df1_320[common_genes_320],df2_320[common_genes_320]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed320, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed320 = np.mean(cos_sim_seed320)
plt.axvline(mean_seed320, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.25, 375, 'mean cosine similarity = '+str(mean_seed320), color='blue', fontsize=10)
plt.text(.25, 325, 'split1 number of genes = '+str(Ngenes_320s1), color='blue', fontsize=10)
plt.text(.25, 275, 'split2 number of genes = '+str(Ngenes_320s2), color='blue', fontsize=10)
plt.text(.25, 225, 'velocities computed with VGENES specified (highly_variable=T)', color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+utv_utvgenes2000, Ngenes='+str(Ngenes_320common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/utvgenes2000_seed320_cos_sim_hist.png')
plt.clf()

# Plot on UMAP
adata_total_res.obs['cos_sim_320'] = cos_sim_seed320
adata_total_res.obs['cos_sim_320'] = pd.DataFrame(adata_total_res.obs['cos_sim_320'])

scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color="cos_sim_320",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/utvgenes2000_seed320_cos_sim_preumap.png")

print("seed320 cosine similarity done!")


