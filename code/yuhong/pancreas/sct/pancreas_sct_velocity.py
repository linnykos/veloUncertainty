### Did not run through all

import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
from torchdiffeq import odeint
from sklearn.metrics.pairwise import cosine_similarity



import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from sctour_misc import *

sct_seed = 615

sct_seed = 615
# https://pytorch.org/docs/stable/notes/randomness.html
torch.manual_seed(sct_seed)
random.seed(sct_seed)
np.random.seed(sct_seed)

""" preprocess code
adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad")
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)  
tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
tnode.train()
adata.obs['ptime'] = tnode.get_time()
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
### 3-tuple of weighted combined latent space, encoder-derived latent space, and ODE-solver-derived latent space
adata.obsm['X_TNODE'] = mix_zs
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
#adata.write(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_sct/pan_sct_preprocess.h5ad")
"""

def train_sct_and_return_tnode_total(adata):
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
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
    return tnode

def train_sct_and_return_tnode(adata):
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
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
    return tnode

# https://stackoverflow.com/questions/28314337/typeerror-sparse-matrix-length-is-ambiguous-use-getnnz-or-shape0-while-usi

###
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_sct/"
total = sc.read(data_folder+"pan_sct_preprocess.h5ad")
total_tnode = train_sct_and_return_tnode_total(total)
total_diff_mat = compute_sctour_velocity(total_tnode, timestep=1/100)
total.layers['velocity'] = total_diff_mat
total.write(data_folder+'res_pancreas_preprocess.h5ad')
total_tnode.save_model(save_dir=data_folder, save_prefix='tnode_pancreas_preprocess')

s1_317 = sc.read(data_folder+"pancreas_seed317_split1_seurat.h5ad")
s1_317_tnode = train_sct_and_return_tnode(s1_317)
s1_317_diff_mat = compute_sctour_velocity(s1_317_tnode, timestep=1/100)

s2_317 = sc.read(data_folder+"pancreas_seed317_split2_seurat.h5ad")
s2_317_tnode = train_sct_and_return_tnode(s2_317)
s2_317_diff_mat = compute_sctour_velocity(s2_317_tnode, timestep=1/100)

# save data
s1_317.layers['velocity'] = s1_317_diff_mat
s1_317.uns['clusters_colors'] = total.uns['clusters_colors'].copy()
s1_317.write(data_folder+'res_pancreas_seed317_split1.h5ad')
s1_317_tnode.save_model(save_dir=data_folder, save_prefix='tnode_pancreas_seed317_split1')
##The first parameter is the directory where you want to save the model, and the second parameter is the prefix for the model name.
# s1_317_tnode.save_model(save_dir='./', save_prefix='pan_sct_tnode_seed317_split1')
# s1_317_tnode = sct.predict.load_model('./pan_sct_tnode_seed317_split1.pth')
# tmp = sct.predict.load_model(data_folder+'tnode_pancreas_seed317_split1.pth')
s2_317.layers['velocity'] = s2_317_diff_mat
s2_317.uns['clusters_colors'] = total.uns['clusters_colors'].copy()
s2_317.write(data_folder+'res_pancreas_seed317_split2.h5ad')
s2_317_tnode.save_model(save_dir=data_folder, save_prefix='tnode_pancreas_seed317_split2')

# cosine similarity histogram
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/"
cos_sim_seed317 = np.diag(cosine_similarity(s1_317_diff_mat,s2_317_diff_mat))
# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black') 
## add mean
plt.axvline(np.mean(cos_sim_seed317), color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-0.2, 500, 'mean cosine similarity = '+str(np.round(np.mean(cos_sim_seed317),4)), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+sct')
plt.savefig(fig_folder+'cos_sim/pan_sct_seed317_hist.png')
plt.clf()


###
s1_320 = sc.read(data_folder+"pancreas_seed320_split1_seurat.h5ad")
s2_320 = sc.read(data_folder+"pancreas_seed320_split2_seurat.h5ad")

s1_320_tnode = train_sct_and_return_tnode(s1_320)
s1_320_diff_mat = compute_sctour_velocity(s1_320_tnode, timestep=1/100)

s2_320_tnode = train_sct_and_return_tnode(s2_320)
s2_320_diff_mat = compute_sctour_velocity(s2_320_tnode, timestep=1/100)

# cosine similarity histogram
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/"
cos_sim_seed320 = np.diag(cosine_similarity(s1_320_diff_mat,s2_320_diff_mat))
# Create histogram
plt.clf()
plt.hist(cos_sim_seed320, bins=30, edgecolor='black') 
## add mean
plt.axvline(np.mean(cos_sim_seed320), color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(-0.2, 500, 'mean cosine similarity = '+str(np.round(np.mean(cos_sim_seed320),4)), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+sct')
plt.savefig(fig_folder+'cos_sim/pan_sct_seed320_hist.png')
plt.clf()



# save data
s1_320.layers['velocity'] = s1_320_diff_mat
s1_320.uns['clusters_colors'] = total.uns['clusters_colors'].copy()
s1_320.write(data_folder+'res_pancreas_seed320_split1.h5ad')
s1_320_tnode.save_model(save_dir=data_folder, save_prefix='tnode_pancreas_seed320_split1')
# tmp = sct.predict.load_model(data_folder+'tnode_pancreas_seed320_split1.pth')
s2_320.layers['velocity'] = s2_320_diff_mat
s2_320.uns['clusters_colors'] = total.uns['clusters_colors'].copy()
s2_320.write(data_folder+'res_pancreas_seed320_split2.h5ad')
s2_320_tnode.save_model(save_dir=data_folder, save_prefix='tnode_pancreas_seed320_split2')



## first check the data object and save it
s1_317.obsm['X_umap'] = total.obsm['X_umap'].copy()
s1_317.obsm['velocity'] = s1_317_diff_mat
sct.vf.plot_vector_field(s1_317, zs_key='X_TNODE', vf_key='velocity', use_rep_neigh='X_TNODE', color='clusters', 
                         show=False, save=fig_folder+'velocity/pan_sct_velocity_seed317_split1.png')

# colors
# Initial arrays
colors_total = total.uns['clusters_colors'] #['#8fbc8f', '#f4a460', '#fdbf6f', '#ff7f00', '#b2df8a', '#1f78b4', '#6a3d9a', '#cab2d6']
clusters_total = np.array(total.obs['clusters'].cat.categories) #['Ductal', 'Ngn3 low EP', 'Ngn3 high EP', 'Pre-endocrine', 'Beta', 'Alpha', 'Delta', 'Epsilon']
new_clusters = np.array(s1_317.obs['clusters'].cat.categories) # ['Alpha', 'Beta', 'Delta', 'Ductal', 'Epsilon', 'Ngn3 high EP', 'Ngn3 low EP', 'Pre-endocrine']
cell_to_color = dict(zip(clusters_total, colors_total))
new_colors = [cell_to_color[cell] for cell in new_clusters]
s1_317.uns['clusters_colors'] = new_colors

## change to scVelo
import scvelo as scv
total = sc.read(data_folder+"pan_sct_preprocess.h5ad")
# umapOriginal
s1_317 = sc.read(data_folder+'res_pancreas_seed317_split1.h5ad')
s1_317.obsm['X_umap'] = total.obsm['X_umap'].copy()
scv.tl.velocity_graph(s1_317)
# WARNING: The neighbor graph has an unexpected format (e.g. computed outside scvelo) 
# or is corrupted (e.g. due to subsetting). Consider recomputing with `pp.neighbors`.
scv.pl.velocity_embedding_stream(s1_317, basis='umap',color="clusters",
                                 save=fig_folder+"velocity/pan_sct_seed317_split1_umapOriginal.png")
# umapCompute
s1_317 = sc.read(data_folder+'res_pancreas_seed317_split1.h5ad')
sc.pp.neighbors(s1_317, use_rep='X_TNODE', n_neighbors=30)
sc.tl.umap(s1_317)
scv.tl.velocity_graph(s1_317)
scv.pl.velocity_embedding_stream(s1_317, basis='umap',color="clusters",
                                 save=fig_folder+"velocity/pan_sct_seed317_split1_umapCompute.png")

def plot_velocity(seed_split):
    data = sc.read(data_folder+'res_pancreas_'+seed_split+'.h5ad')
    data.uns['clusters_colors'] = new_colors # will delete this line
    data.obsm['X_umap'] = total.obsm['X_umap'].copy()
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis='umap',color="clusters",save=fig_folder+"velocity/pan_sct_"+seed_split+"_umapOriginal.png")
    data = sc.read(data_folder+'res_pancreas_'+seed_split+'.h5ad')
    data.uns['clusters_colors'] = new_colors # will delete this line
    sc.pp.neighbors(data, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(data)
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis='umap',color="clusters",save=fig_folder+"velocity/pan_sct_"+seed_split+"_umapCompute.png")

plot_velocity("seed317_split1")
plot_velocity("seed317_split2")
plot_velocity("seed320_split1")
plot_velocity("seed320_split2")


s1_317 = sc.read(data_folder+'res_pancreas_seed317_split1.h5ad')
s2_317 = sc.read(data_folder+'res_pancreas_seed317_split2.h5ad')
s1_320 = sc.read(data_folder+'res_pancreas_seed320_split1.h5ad')
s2_320 = sc.read(data_folder+'res_pancreas_seed320_split2.h5ad')


total = sc.read(data_folder+'res_pancreas_preprocess.h5ad')
sc.pp.neighbors(total, use_rep='X_TNODE', n_neighbors=30)
# do this because of the error:
# ValueError: Your neighbor graph seems to be corrupted. Consider recomputing via pp.neighbors.
scv.tl.velocity_graph(total)
total.obs['cos_sim_317'] = np.diag(cosine_similarity(s1_317.layers['velocity'],s2_317.layers['velocity']))
total.obs['cos_sim_317'] = pd.DataFrame(total.obs['cos_sim_317'])
scv.pl.velocity_embedding_stream(total, basis='umap',color="cos_sim_317",cmap='coolwarm',
                                 save=fig_folder+"cos_sim/pan_sct_seed317_cos_similarity.png")
total.obs['cos_sim_320'] = np.diag(cosine_similarity(s1_320.layers['velocity'],s2_320.layers['velocity']))
total.obs['cos_sim_320'] = pd.DataFrame(total.obs['cos_sim_320'])
scv.pl.velocity_embedding_stream(total, basis='umap',color="cos_sim_320",cmap='coolwarm',
                                 save=fig_folder+"cos_sim/pan_sct_seed320_cos_similarity.png")

scv.pl.velocity_embedding_stream(total, basis='umap',color="clusters", save=fig_folder+"velocity/pan_sct_preprocess.png")


