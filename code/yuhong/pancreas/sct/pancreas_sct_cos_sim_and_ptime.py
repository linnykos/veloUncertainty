import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
import torch
import random

sct_seed = 615
# https://pytorch.org/docs/stable/notes/randomness.html

def train_sct(adata):
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
    return adata

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_sct/"
total = sc.read(data_folder+"pan_sct_preprocess.h5ad")

s1_317 = sc.read(data_folder+"pancreas_seed317_split1_seurat.h5ad")
s1_317_res = train_sct(s1_317)
s1_317_res.write(filename=data_folder+"pancreas_seed317_split1.h5ad")

s2_317 = sc.read(data_folder+"pancreas_seed317_split2_seurat.h5ad")
s2_317_res = train_sct(s2_317)
s2_317_res.write(filename=data_folder+"pancreas_seed317_split2.h5ad")

s1_320 = sc.read(data_folder+"pancreas_seed320_split1_seurat.h5ad")
s1_320_res = train_sct(s1_320)
s1_320_res.write(filename=data_folder+"pancreas_seed320_split1.h5ad")

s2_320 = sc.read(data_folder+"pancreas_seed320_split2_seurat.h5ad")
s2_320_res = train_sct(s2_320)
s2_320_res.write(filename=data_folder+"pancreas_seed320_split2.h5ad")

s1_317_res.uns['clusters_colors']=total.uns['clusters_colors']
s2_317_res.uns['clusters_colors']=total.uns['clusters_colors']
s1_320_res.uns['clusters_colors']=total.uns['clusters_colors']
s2_320_res.uns['clusters_colors']=total.uns['clusters_colors']

#########################################################
### visualization
fig_umap_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/umap/"
def plot_vf_umap(res,name):
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(res, color='clusters', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    sc.pl.umap(res, color='ptime', ax=axs[1], show=False, frameon=False)
    sct.vf.plot_vector_field(res, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='clusters', 
                            show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                            save=fig_umap_folder+name)
s1_317_res.obsm['X_umap'] = total.obsm['X_umap'].copy()
s2_317_res.obsm['X_umap'] = total.obsm['X_umap'].copy()
s1_320_res.obsm['X_umap'] = total.obsm['X_umap'].copy()
s2_320_res.obsm['X_umap'] = total.obsm['X_umap'].copy()

plot_vf_umap(res=s1_317_res, name="pancreas_sct_vf_seed317_split1_preumap.png")
plot_vf_umap(res=s2_317_res, name="pancreas_sct_vf_seed317_split2_preumap.png")
plot_vf_umap(res=s1_320_res, name="pancreas_sct_vf_seed320_split1_preumap.png")
plot_vf_umap(res=s2_320_res, name="pancreas_sct_vf_seed320_split2_preumap.png")

def compute_umap_and_plot_vf(res,name):
    del res.obsm['X_umap']
    res = res[np.argsort(res.obs['ptime'].values), :]
    sc.pp.neighbors(res, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(res)
    plot_vf_umap(res=res, name=name)

compute_umap_and_plot_vf(res=s1_317_res,name="pancreas_sct_vf_seed317_split1_computeumap.png")
compute_umap_and_plot_vf(res=s2_317_res, name="pancreas_sct_vf_seed317_split2_computeumap.png")
compute_umap_and_plot_vf(res=s1_320_res, name="pancreas_sct_vf_seed320_split1_computeumap.png")
compute_umap_and_plot_vf(res=s2_320_res, name="pancreas_sct_vf_seed320_split2_computeumap.png")


#########################################################
### cosine similarity
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_sct/"
total = sc.read(data_folder+"pan_sct_preprocess.h5ad")
s1_317_res = sc.read(data_folder+"pancreas_seed317_split1.h5ad")
s2_317_res = sc.read(data_folder+"pancreas_seed317_split2.h5ad")
s1_320_res = sc.read(data_folder+"pancreas_seed320_split1.h5ad")
s2_320_res = sc.read(data_folder+"pancreas_seed320_split2.h5ad")

cos_sim_seed317 = np.diag(cosine_similarity(s1_317_res.obsm['X_VF'],s2_317_res.obsm['X_VF']))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black') 
## add mean
mean_seed317 = np.mean(cos_sim_seed317)
plt.axvline(mean_seed317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.3, 500, 'mean cosine similarity = '+str(np.round(mean_seed317,4)), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+sct')
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/pan_sct_seed317_cos_similarity_hist.png')
plt.clf()

total.obs['cos_sim_317'] = cos_sim_seed317
total.obs['cos_sim_317'] = pd.DataFrame(total.obs['cos_sim_317'])

sct.vf.plot_vector_field(total, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='cos_sim_317', cmap='coolwarm',
                         legend_loc='none', frameon=False, size=100, alpha=0.2, 
                         save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/pan_sct_seed317_cos_sim.png")


## seed=320
cos_sim_seed320 = np.diag(cosine_similarity(s1_320_res.obsm['X_VF'],s2_320_res.obsm['X_VF']))

plt.clf()
plt.hist(cos_sim_seed320, bins=30, edgecolor='black') 
## add mean
mean_seed320 = np.mean(cos_sim_seed320)
plt.axvline(mean_seed320, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.9, 800, 'mean cosine similarity = '+str(np.round(mean_seed320,4)), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed320)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+sct')
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/pan_sct_seed320_cos_similarity_hist.png')
plt.clf()

total.obs['cos_sim_320'] = cos_sim_seed320
total.obs['cos_sim_320'] = pd.DataFrame(total.obs['cos_sim_320'])

sct.vf.plot_vector_field(total, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='cos_sim_320', cmap='coolwarm',
                         legend_loc='none', frameon=False, size=100, alpha=0.2, 
                         save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/pan_sct_seed320_cos_sim.png")


#########################################################
### pseudotime
## pseudotime in these objects have been reordered for the purpose of plotting
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_sct/"
total = sc.read(data_folder+"pan_sct_preprocess.h5ad")
s1_317_res = sc.read(data_folder+"pancreas_seed317_split1.h5ad")
s2_317_res = sc.read(data_folder+"pancreas_seed317_split2.h5ad")
s1_320_res = sc.read(data_folder+"pancreas_seed320_split1.h5ad")
s2_320_res = sc.read(data_folder+"pancreas_seed320_split2.h5ad")

np.corrcoef(s1_317_res.obs['ptime'],s2_317_res.obs['ptime']) 
np.corrcoef(s1_320_res.obs['ptime'],s2_320_res.obs['ptime']) 

np.corrcoef(s1_317_res.obs['ptime'],total.obs['ptime']) 
np.corrcoef(total.obs['ptime'],s2_317_res.obs['ptime']) 
np.corrcoef(s1_320_res.obs['ptime'],total.obs['ptime']) 
np.corrcoef(total.obs['ptime'],s2_320_res.obs['ptime']) 

out_fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/sctour/Jun13_pseudotime/"
cell_types = total.obs['clusters']
colors = dict(zip(['Ductal', 'Ngn3 low EP', 'Ngn3 high EP', 'Pre-endocrine', 'Beta','Alpha', 'Delta', 'Epsilon'],
                  ['#8fbc8f', '#f4a460', '#fdbf6f', '#ff7f00', '#b2df8a', '#1f78b4','#6a3d9a', '#cab2d6']))
def ptime_scatter_plot(s1,s2,seed,method,name,xlab,ylab,data="pan"):
    df = pd.DataFrame({'split1':np.array(s1.obs['ptime']), 'split2':np.array(s2.obs['ptime']),
                       'cell_types':np.array(cell_types)})
    corr = np.round(np.corrcoef(s1.obs['ptime'],s2.obs['ptime']),3)[0,1]
    print(corr)
    for category, color in colors.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors))
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime '+name+' (seed='+str(seed)+', corr='+str(corr)+')')
    plt.savefig(out_fig_folder+data+"_"+method+"_seed"+str(seed)+"_pseudotime"+name+".png")
    plt.close()

ptime_scatter_plot(s1=s1_317_res,s2=s2_317_res,seed=317,method="sct",name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=s1_317_res,s2=total,seed=317,method="sct",name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=s2_317_res,s2=total,seed=317,method="sct",name="split2vstotal",xlab="split2",ylab="total")

ptime_scatter_plot(s1=s1_320_res,s2=s2_320_res,seed=320,method="sct",name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=s1_320_res,s2=total,seed=320,method="sct",name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=s2_320_res,s2=total,seed=320,method="sct",name="split2vstotal",xlab="split2",ylab="total")


