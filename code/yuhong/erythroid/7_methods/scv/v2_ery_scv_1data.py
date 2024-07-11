import datetime
import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/scv/" 

print_message_with_time("########### Start to read total, split1 and split2")
total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
adata_split1 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split1_allgenes.h5ad') # 9815 × 53801
adata_split2 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split2_allgenes.h5ad')
gene_names = total.var.index.copy()
S_mat_split1 = adata_split1.layers['spliced'].copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()

def scv_compute_velocity(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    ### batch correction
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("Batch correction done!")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

scv_compute_velocity(total) # running this
scv_compute_velocity(adata_split1)
scv_compute_velocity(adata_split2)

positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions_split1 = [positions_dict[gene] for gene in adata_split1.var.index]
positions_split2 = [positions_dict[gene] for gene in adata_split2.var.index]

S_mat_split1 = S_mat_split1[:,positions_split1]
U_mat_split1 = U_mat_split1[:,positions_split1]
adata_split1.layers['spliced_original'] = S_mat_split1[:,positions_split1] # looks like i did not do this actually
adata_split1.layers['unspliced_original'] = U_mat_split1[:,positions_split1]
adata_split2.layers['spliced_original'] = S_mat_split2[:,positions_split2]
adata_split2.layers['unspliced_original'] = U_mat_split2[:,positions_split2]


import numpy as np
common_genes_filter = np.intersect1d(np.array(adata_split1.var.index), np.array(adata_split2.var.index))
common_genes_filter.shape # 1481 common genes
np.intersect1d(np.array(adata_split1.var.index), np.array(total.var.index)).shape # 589 -> 1335
np.intersect1d(np.array(adata_split2.var.index), np.array(total.var.index)).shape # 580 -> 1328

#adata_split1.layers["velocity_rmNA"] = np.nan_to_num(adata_split1.layers['velocity'], nan=0)
#adata_split2.layers["velocity_rmNA"] = np.nan_to_num(adata_split2.layers['velocity'], nan=0)
np.sum(~np.isnan(adata_split1.layers['velocity'][0])) # Ngenes in split1 for velocity computation=9 -> 311
np.sum(~np.isnan(adata_split2.layers['velocity'][0])) # Ngenes in split2 for velocity computation=9 -> 315

velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
common_genes_velocity.shape # =1 -> 209

# write data
total.write_h5ad(data_folder+'v2_erythroid/scv/adata_ery_scv_total_v2.h5ad')
adata_split1.write_h5ad(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split1_v2.h5ad')
adata_split2.write_h5ad(data_folder+'v2_erythroid/scv/adata_ery_scv_seed317_split2_v2.h5ad')

######################################################
## plot velocities
raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
scv.pl.velocity_embedding_stream(total, basis='umap',color="celltype",save=fig_folder+"velocity/ery_scv_total_umapCompute.png")
scv.pl.velocity_embedding_stream(total, basis='umap',color="sequencing.batch",save=fig_folder+"velocity/ery_scv_total_umapCompute_byBatch.png")

def plot_velocities_scv(split, adata_total, fig_info):
    adata=split.copy()
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="celltype",save=fig_folder+"velocity/ery_scv_"+fig_info+"_umapCompute.png")
    adata.obsm['X_umap'] = adata_total.obsm['X_umap'].copy()
    scv.pl.velocity_embedding_stream(adata, basis='umap',color="celltype",save=fig_folder+"velocity/ery_scv_"+fig_info+"_umapOriginal.png")

plot_velocities_scv(total,raw,"total")
plot_velocities_scv(adata_split1,raw,"seed317_split1")
plot_velocities_scv(adata_split2,raw,"seed317_split2")


import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
velo_df1 = pd.DataFrame(adata_split1.layers['velocity'], columns=adata_split1.var.index.tolist())
velo_df2 = pd.DataFrame(adata_split2.layers['velocity'], columns=adata_split2.var.index.tolist())
cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))
np.quantile(cos_sim,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

# will be different
#cs1=np.diag(cosine_similarity(np.nan_to_num(adata_split1.layers['velocity'], nan=0),np.nan_to_num(adata_split2.layers['velocity'], nan=0)))
#np.quantile(cs1,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])

import matplotlib.pyplot as plt
# cosine similarity
def plot_cosine_similarity_scv(cos_sim, total, seed_split=317,text_x=None,text_y=None):
    # histogram
    plt.clf()
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='dimgray',color='powderblue') 
    max_frequency = np.max(counts)
    if text_x is None:
        text_x = np.quantile(cos_sim,[.05])[0]
    if text_y is None:
        text_y = max_frequency/2
    plt.axvline(np.mean(cos_sim), color='salmon', linestyle='dashed', linewidth=1.5) ## add mean
    plt.text(text_x, text_y, 'mean = '+str(np.round(np.mean(cos_sim),5)), color='navy', fontsize=11)
    plt.xlabel('cosine similarity (seed'+str(seed_split)+')')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, ery+scv')
    plt.savefig(fig_folder+'cos_sim/seed'+str(seed_split)+'_cos_sim_hist.png')
    plt.clf()
    # umap
    total.obs['cos_sim'] = cos_sim
    scv.pl.velocity_embedding_stream(total, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                 save=fig_folder+"cos_sim/seed"+str(seed_split)+"_cos_sim_umapCompute.png")

plot_cosine_similarity_scv(cos_sim,total)

scv.tl.velocity_confidence(total)
scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[1, 100],
               save=fig_folder+"velo_conf/total_veloConf_umapCompute.png")
scv.pl.scatter(total, color='cos_sim', cmap='coolwarm', perc=[1, 100],
               save=fig_folder+"cos_sim/seed317_cos_sim_umapCompute_scatter.png")


