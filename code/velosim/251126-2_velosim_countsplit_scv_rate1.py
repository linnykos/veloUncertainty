# velosim
import pandas as pd
import numpy as np
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/velosim/251126-1_rate1/'
method = 'scv'
dataset_short = 'velosim1126_rate1'
dataset_long = 'velosim1126_rate1'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/velosim/251126-1_rate1/'+method+'/'

S = pd.read_csv(data_folder+'velosim_counts_s.csv')
U = pd.read_csv(data_folder+'velosim_counts_u.csv')
overdisp_S = np.array(pd.read_csv(data_folder+'velosim_overdisp_S.csv')['x'])
overdisp_U = np.array(pd.read_csv(data_folder+'velosim_overdisp_U.csv')['x'])

split_seed = 317

np.random.seed(split_seed)
s1, s2  = countsplit(S,overdisps=overdisp_S)
u1, u2  = countsplit(U,overdisps=overdisp_U)

velo_true = pd.read_csv(data_folder+'velosim_velocity.csv')
ptime = pd.read_csv(data_folder+'velosim_pseudo_time.csv')

# create adata objects
import anndata as ad
import scipy.sparse as sp

def create_adata(S, U, velo_true):
    adata = ad.AnnData(X = sp.csr_matrix(S))
    adata.layers["spliced"] = sp.csr_matrix(S)
    adata.layers["unspliced"] = sp.csr_matrix(U)
    adata.layers['velocity_true'] = velo_true
    return adata

# total
adata = create_adata(S, U, velo_true)
# split1
adata1 = create_adata(s1, u1, velo_true)
# split2
adata2 = create_adata(s2, u2, velo_true)


def run_scv_velosim(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata,n_jobs=8)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata,n_jobs=8)

run_scv_velosim(adata1)
run_scv_velosim(adata2)
run_scv_velosim(adata)

# make plots
import scvelo as scv
import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

def make_velosim_plots(fig_folder, adata, adata1, adata2, dataset_short, method, split_seed):
    adata.obsm['X_umapOriginal'] = adata.obsm['X_umap'].copy()
    adata1.obsm['X_umapOriginal'] = adata1.obsm['X_umap'].copy()
    adata2.obsm['X_umapOriginal'] = adata2.obsm['X_umap'].copy()
    plot_velocity(adata_in=adata, fig_folder=fig_folder, data_version="total", dataset=dataset_short, method=method, split_seed=split_seed)
    plot_velocity(adata_in=adata1, fig_folder=fig_folder, data_version="split1", dataset=dataset_short, method=method, split_seed=split_seed)
    plot_velocity(adata_in=adata2, fig_folder=fig_folder, data_version="split2", dataset=dataset_short, method=method, split_seed=split_seed)
    plot_cosine_similarity(adata_split1=adata1,adata_split2=adata2,adata_total=adata,
                           dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    plot_cosine_similarity_withRef(adata1,adata2,adata,dataset_short,method,fig_folder,split_seed)
    # median cos_sim=-0.2162, mean cos_sim=-0.2196
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_confidence(adata1)
    scv.tl.velocity_confidence(adata1)
    plot_veloConf_and_cosSim(adata,adata1,adata2,dataset_short,method,fig_folder,split_seed)
    plot_veloConf_hist(adata,dataset_short,method,fig_folder,split_seed)
    if method=='scv':
        mask = ~( np.all(np.isnan(adata.layers['velocity']), axis=0) )
        vel_filtered = adata.layers['velocity'][:, mask]
        vel_true_filtered = adata.layers['velocity_true'][:, mask]
        cos_sim = np.diag(cosine_similarity(vel_filtered, vel_true_filtered))
        print('mean cos_sim with true: '+str(np.mean( cos_sim ) ) )
        print('median cos_sim with true: ', str(np.median( cos_sim ) ) )
    else:
        cos_sim = np.diag(cosine_similarity(adata.layers['velocity'],adata.layers['velocity_true']))
        print('mean cos_sim with true: '+str(np.mean( cos_sim ) ) )
        print('median cos_sim with true: ', str(np.median( cos_sim ) ) )
    # plot 
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='gainsboro',color='powderblue') 
    max_frequency = np.max(counts)
    text_x = np.quantile(cos_sim,[.05])[0]
    text_y = max_frequency/5
    plt.axvline(np.mean(cos_sim), color='brown', linestyle='dashed', linewidth=1.5) ## add mean
    plt.axvline(np.median(cos_sim), color='peru', linestyle='dashed', linewidth=1.5) ## add median
    plt.text(text_x,text_y*2.5,'mean='+str(np.round(np.mean(cos_sim),4)), color='firebrick', fontsize=11)
    plt.text(text_x,text_y*3,'median='+str(np.round(np.median(cos_sim),4)), color='sienna', fontsize=11)
    plt.xlabel('cosine similarity')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity (computed vs true), '+dataset_short+'+'+method)
    plt.savefig(fig_folder+'metric/'+dataset_short+'+'+method+'_cos_sim_hist_comp_vs_true.png')
    plt.clf()

make_velosim_plots(fig_folder, adata, adata1, adata2, dataset_short, method, split_seed)





