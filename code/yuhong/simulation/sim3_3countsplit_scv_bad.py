import scvelo as scv
import scanpy as sc
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import matplotlib.pyplot as plt

data = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3bad.h5ad")
scv.pp.log1p(data)
sc.pp.pca(data)
sc.pp.neighbors(data, n_pcs=30, n_neighbors=30)
scv.pp.moments(data, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(data,t_max=25)
scv.tl.velocity_graph(data,vkey="true_velocity")

# seed317, split1
s1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3bad_seed317_split1_seurat.h5ad")
s1.obs['true_t'] = data.obs['true_t'] 
scv.pp.log1p(s1)
sc.pp.pca(s1)
sc.pp.neighbors(s1, n_pcs=30, n_neighbors=30)
scv.pp.moments(s1, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(s1)
scv.tl.velocity(s1, mode="dynamical")
scv.tl.velocity_graph(s1)
scv.pl.velocity_embedding_stream(s1,color="true_t", basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_317split1_vcompute.png")

# seed317, split2
s2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3bad_seed317_split2_seurat.h5ad")
s2.obs['true_t'] = data.obs['true_t'] 
scv.pp.log1p(s2)
sc.pp.pca(s2)
sc.pp.neighbors(s2, n_pcs=30, n_neighbors=30)
scv.pp.moments(s2, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(s2)
scv.tl.velocity(s2, mode="dynamical")
scv.tl.velocity_graph(s2)
scv.pl.velocity_embedding_stream(s2,color="true_t", basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_317split2_vcompute.png")

# seed317, cosine similarity
s1.layers["velocity_rmNA"] = np.nan_to_num(s1.layers['velocity'], nan=0)
s2.layers["velocity_rmNA"] = np.nan_to_num(s2.layers['velocity'], nan=0)
cos_sim = np.diag(cosine_similarity(s1.layers["velocity_rmNA"], s2.layers["velocity_rmNA"]))
# Create histogram
plt.clf()
plt.hist(cos_sim, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_317 = np.mean(cos_sim)
Ngenes_317s1 = np.sum(~np.isnan(s1.layers['velocity'][0]))
Ngenes_317s2 = np.sum(~np.isnan(s2.layers['velocity'][0]))
Ngenes_317common = np.sum(np.isnan(s1.layers["velocity"][0] + s2.layers["velocity"][0])==0)
plt.axvline(mean_317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.1, 80, 'mean cosine similarity = '+str(mean_317), color='blue', fontsize=10)
plt.text(.1, 70, 'split1 number of genes'+str(Ngenes_317s1), color='blue', fontsize=10)
plt.text(.1, 60, 'split2 number of genes'+str(Ngenes_317s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, sim3bad, Ngenes='+str(Ngenes_317common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_cos_sim_317_hist.png')
plt.clf()

data.obs["cos_sim"] = cos_sim
data.obs["cos_sim"] = pd.DataFrame(data.obs["cos_sim"])
scv.pl.velocity_embedding_stream(data, vkey="true_velocity",basis='pca',color="cos_sim",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_cos_sim_317_pca.png")

###
# seed320, split1
s1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3bad_seed320_split1_seurat.h5ad")
s1.obs['true_t'] = data.obs['true_t'] 
scv.pp.log1p(s1)
sc.pp.pca(s1)
sc.pp.neighbors(s1, n_pcs=30, n_neighbors=30)
scv.pp.moments(s1, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(s1)
scv.tl.velocity(s1, mode="dynamical")
scv.tl.velocity_graph(s1)
scv.pl.velocity_embedding_stream(s1,color="true_t", basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_320split1_vcompute.png")

# seed320, split2
s2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3bad_seed320_split2_seurat.h5ad")
s2.obs['true_t'] = data.obs['true_t'] 
scv.pp.log1p(s2)
sc.pp.pca(s2)
sc.pp.neighbors(s2, n_pcs=30, n_neighbors=30)
scv.pp.moments(s2, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(s2)
scv.tl.velocity(s2, mode="dynamical")
scv.tl.velocity_graph(s2)
scv.pl.velocity_embedding_stream(s2,color="true_t", basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_320split2_vcompute.png")

# seed320, cosine similarity
s1.layers["velocity_rmNA"] = np.nan_to_num(s1.layers['velocity'], nan=0)
s2.layers["velocity_rmNA"] = np.nan_to_num(s2.layers['velocity'], nan=0)
cos_sim = np.diag(cosine_similarity(s1.layers["velocity_rmNA"], s2.layers["velocity_rmNA"]))
# Create histogram
plt.clf()
plt.hist(cos_sim, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_320 = np.mean(cos_sim)
Ngenes_320s1 = np.sum(~np.isnan(s1.layers['velocity'][0]))
Ngenes_320s2 = np.sum(~np.isnan(s2.layers['velocity'][0]))
Ngenes_320common = np.sum(np.isnan(s1.layers["velocity"][0] + s2.layers["velocity"][0])==0)
plt.axvline(mean_317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.1, 80, 'mean cosine similarity = '+str(mean_320), color='blue', fontsize=10)
plt.text(.1, 70, 'split1 number of genes'+str(Ngenes_320s1), color='blue', fontsize=10)
plt.text(.1, 60, 'split2 number of genes'+str(Ngenes_320s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, sim3bad, Ngenes='+str(Ngenes_320common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_cos_sim_320_hist.png')
plt.clf()

data.obs["cos_sim"] = cos_sim
data.obs["cos_sim"] = pd.DataFrame(data.obs["cos_sim"])
scv.pl.velocity_embedding_stream(data, vkey="true_velocity",basis='pca',color="cos_sim",cmap='coolwarm',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_cos_sim_320_pca.png")


