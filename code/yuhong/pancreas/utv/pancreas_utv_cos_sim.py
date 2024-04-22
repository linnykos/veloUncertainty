##### cosine similarity
import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata_total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_total.h5ad')
split1_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed317_split1.h5ad')
split2_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed317_split2.h5ad')

#split1_seed317_res.layers["velocity_rmNA"] = np.nan_to_num(split1_seed317_res.layers['velocity'], nan=0)
#split2_seed317_res.layers["velocity_rmNA"] = np.nan_to_num(split2_seed317_res.layers['velocity'], nan=0)
common_genes_317 = np.intersect1d(split1_seed317_res.var['features'], split2_seed317_res.var['features'])
#indices_seed317 = split1_seed317_res.var['features'][split1_seed317_res.var['features'].isin(set(common_genes_317))].index.tolist()
Ngenes_317s1 = len(split1_seed317_res.var['features'])
Ngenes_317s2 = len(split2_seed317_res.var['features'])
Ngenes_317common = len(common_genes_317)

# Here assuming order of genes in layers['velocity'] is the same as in var['features']
### split1_seed317_res.layers['velocity'].shape = (3696, 734)
### split2_seed317_res.layers['velocity'].shape = (3696, 731)
df1 = pd.DataFrame(split1_seed317_res.layers['velocity'], columns=split1_seed317_res.var['features'].tolist())
df2 = pd.DataFrame(split2_seed317_res.layers['velocity'], columns=split2_seed317_res.var['features'].tolist())
cos_sim_seed317 = np.diag(cosine_similarity(df1[common_genes_317],df2[common_genes_317]))
#cos_sim_seed317 = np.diag(cosine_similarity(df1[indices_seed317],df2[indices_seed317]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed317 = np.mean(cos_sim_seed317)
plt.axvline(mean_seed317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-.75, 320, 'mean cosine similarity = '+str(mean_seed317), color='blue', fontsize=10)
plt.text(-.75, 300, 'split1 number of genes = '+str(Ngenes_317s1), color='blue', fontsize=10)
plt.text(-.75, 280, 'split2 number of genes = '+str(Ngenes_317s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+utv, Ngenes='+str(Ngenes_317common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_cos_similarity_hist.png')
plt.clf()

## Plot on UMAP
adata_total_res.obs['cos_sim_317'] = cos_sim_seed317
adata_total_res.obs['cos_sim_317'] = pd.DataFrame(adata_total_res.obs['cos_sim_317'])
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_cos_similarity.png")
del adata_total_res.obs['cos_sim_317']
print("seed317 cosine similarity done!")

#############################
split1_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed320_split1.h5ad')
split2_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed320_split2.h5ad')

Ngenes_320s1 = len(split1_seed320_res.var['features']) # 741
Ngenes_320s2 = len(split2_seed320_res.var['features']) # 744
common_genes_320 = np.intersect1d(split1_seed320_res.var['features'], split2_seed320_res.var['features'])
Ngenes_320common = len(common_genes_320) # 654
# indices_seed320 = split1_seed320_res.var['features'][split1_seed320_res.var['features'].isin(set(common_genes_320))].index.tolist()

df1 = pd.DataFrame(split1_seed320_res.layers['velocity'], columns=split1_seed320_res.var['features'].tolist())
df2 = pd.DataFrame(split2_seed320_res.layers['velocity'], columns=split2_seed320_res.var['features'].tolist())
cos_sim_seed320 = np.diag(cosine_similarity(df1[common_genes_320],df2[common_genes_320]))
#cos_sim_seed320 = np.diag(cosine_similarity(df1[indices_seed320],df2[indices_seed320]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed320, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed320 = np.mean(cos_sim_seed320)
plt.axvline(mean_seed320, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-.75, 320, 'mean cosine similarity = '+str(mean_seed320), color='blue', fontsize=10)
plt.text(-.75, 300, 'split1 number of genes = '+str(Ngenes_320s1), color='blue', fontsize=10)
plt.text(-.75, 280, 'split2 number of genes = '+str(Ngenes_320s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+utv, Ngenes='+str(Ngenes_320common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed320_cos_similarity_hist.png')
plt.clf()

# Plot on UMAP
adata_total_res.obs['cos_sim_320'] = cos_sim_seed320
adata_total_res.obs['cos_sim_320'] = pd.DataFrame(adata_total_res.obs['cos_sim_320'])

scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color="cos_sim_320",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed320_cos_similarity.png")
print("seed320 cosine similarity done!")


