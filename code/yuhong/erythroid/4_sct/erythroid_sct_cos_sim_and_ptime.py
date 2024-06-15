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

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_sct/"
total = sc.read(data_folder+"ery_sct_preprocess.h5ad")

s1_317 = sc.read(data_folder+"erythroid_seed317_split1_seurat.h5ad")
s1_317_res = train_sct(s1_317)
s1_317_res.write(filename=data_folder+"erythroid_seed317_split1.h5ad")

s2_317 = sc.read(data_folder+"erythroid_seed317_split2_seurat.h5ad")
s2_317_res = train_sct(s2_317)
s2_317_res.write(filename=data_folder+"erythroid_seed317_split2.h5ad")

s1_320 = sc.read(data_folder+"erythroid_seed320_split1_seurat.h5ad")
s1_320_res = train_sct(s1_320)
s1_320_res.write(filename=data_folder+"erythroid_seed320_split1.h5ad")

s2_320 = sc.read(data_folder+"erythroid_seed320_split2_seurat.h5ad")
s2_320_res = train_sct(s2_320)
s2_320_res.write(filename=data_folder+"erythroid_seed320_split2.h5ad")

#########################################
### pseudotime

np.corrcoef(s1_317_res.obs['ptime'],s2_317_res.obs['ptime']) 
np.corrcoef(s1_320_res.obs['ptime'],s2_320_res.obs['ptime']) 

np.corrcoef(s1_317_res.obs['ptime'],total.obs['ptime']) 
np.corrcoef(total.obs['ptime'],s2_317_res.obs['ptime']) 
np.corrcoef(s1_320_res.obs['ptime'],total.obs['ptime']) 
np.corrcoef(total.obs['ptime'],s2_320_res.obs['ptime']) 

out_fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/sctour/pseudotime/"
cell_types = total.obs['celltype']
colors = dict(zip(['Blood progenitors 1', 'Blood progenitors 2', 'Erythroid1', 'Erythroid2', 'Erythroid3'],
                  ['#f9decf', '#c9a997', '#C72228', '#f79083', '#EF4E22']))
def ptime_scatter_plot(s1,s2,seed,method,name,xlab,ylab,data="ery"):
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


#########################################
### cosine similarity
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_sct/"
total = sc.read(data_folder+"ery_sct_preprocess.h5ad")
s1_317_res = sc.read(data_folder+"erythroid_seed317_split1.h5ad")
s2_317_res = sc.read(data_folder+"erythroid_seed317_split2.h5ad")
s1_320_res = sc.read(data_folder+"erythroid_seed320_split1.h5ad")
s2_320_res = sc.read(data_folder+"erythroid_seed320_split2.h5ad")

cos_sim_seed317 = np.diag(cosine_similarity(s1_317_res.obsm['X_VF'],s2_317_res.obsm['X_VF']))
cos_sim_seed320 = np.diag(cosine_similarity(s1_320_res.obsm['X_VF'],s2_320_res.obsm['X_VF']))

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/sctour/cos_sim/"
def plot_cos_sim(cos_sim,text_x,text_y,fig_name_hist,fig_name_umap,seed):
    plt.clf()
    plt.hist(cos_sim, bins=30, edgecolor='black') 
    mean_cos_sim = np.mean(cos_sim)  ## add mean
    plt.axvline(mean_cos_sim, color='red', linestyle='dashed', linewidth=1)
    plt.text(text_x, text_y, 'mean cosine similarity = '+str(np.round(mean_cos_sim,4)), color='blue', fontsize=10)
    ## add labels and title
    plt.xlabel('cosine similarity (seed='+str(seed)+')')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, ery+sct')
    plt.savefig(fig_folder+fig_name_hist)
    plt.clf()
    total.obs['cos_sim_'+str(seed)] = cos_sim
    total.obs['cos_sim_'+str(seed)] = pd.DataFrame(total.obs['cos_sim_'+str(seed)])
    sct.vf.plot_vector_field(total, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='cos_sim_'+str(seed), cmap='coolwarm',
                         legend_loc='none', frameon=False, size=100, alpha=0.2, save=fig_folder+fig_name_umap)


plot_cos_sim(cos_sim=cos_sim_seed317,text_x=.4,text_y=600,fig_name_hist='ery_sct_seed317_cos_sim_hist.png',
             fig_name_umap='ery_sct_seed317_cos_sim.png',seed=317)
plot_cos_sim(cos_sim=cos_sim_seed320,text_x=0,text_y=200,fig_name_hist='ery_sct_seed320_cos_sim_hist.png',
             fig_name_umap='ery_sct_seed320_cos_sim.png',seed=320)
## RuntimeWarning: invalid value encountered in divide
### cos_sim = np.einsum("ij, j", dZ, V[i]) / (l2_norm(dZ, axis = 1) * l2_norm(V[i]))

#########################################################
### UMAP visualization
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_sct/"
total = sc.read(data_folder+"ery_sct_preprocess.h5ad")
s1_317_res = sc.read(data_folder+"erythroid_seed317_split1.h5ad")
s2_317_res = sc.read(data_folder+"erythroid_seed317_split2.h5ad")
s1_320_res = sc.read(data_folder+"erythroid_seed320_split1.h5ad")
s2_320_res = sc.read(data_folder+"erythroid_seed320_split2.h5ad")

fig_umap_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/sctour/umap/"
def plot_vf_umap(res,name):
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(res, color='celltype', ax=axs[0], legend_loc='on data', show=False, frameon=False)
    sc.pl.umap(res, color='ptime', ax=axs[1], show=False, frameon=False)
    sct.vf.plot_vector_field(res, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='celltype', 
                            show=False, ax=axs[2], legend_loc='none', frameon=False, size=100, alpha=0.2, 
                            save=fig_umap_folder+name)
s1_317_res.obsm['X_umap'] = total.obsm['X_umap'].copy()
s2_317_res.obsm['X_umap'] = total.obsm['X_umap'].copy()
s1_320_res.obsm['X_umap'] = total.obsm['X_umap'].copy()
s2_320_res.obsm['X_umap'] = total.obsm['X_umap'].copy()

plot_vf_umap(res=s1_317_res, name="ery_sct_vf_seed317_split1_preumap.png")
plot_vf_umap(res=s2_317_res, name="ery_sct_vf_seed317_split2_preumap.png")
plot_vf_umap(res=s1_320_res, name="ery_sct_vf_seed320_split1_preumap.png")
plot_vf_umap(res=s2_320_res, name="ery_sct_vf_seed320_split2_preumap.png")

def compute_umap_and_plot_vf(res,name):
    del res.obsm['X_umap']
    res = res[np.argsort(res.obs['ptime'].values), :]
    sc.pp.neighbors(res, use_rep='X_TNODE', n_neighbors=30)
    sc.tl.umap(res)
    plot_vf_umap(res=res, name=name)

compute_umap_and_plot_vf(res=s1_317_res,name="ery_sct_vf_seed317_split1_computeumap.png")
compute_umap_and_plot_vf(res=s2_317_res, name="ery_sct_vf_seed317_split2_computeumap.png")
compute_umap_and_plot_vf(res=s1_320_res, name="ery_sct_vf_seed320_split1_computeumap.png")
compute_umap_and_plot_vf(res=s2_320_res, name="ery_sct_vf_seed320_split2_computeumap.png")



