import scvelo as scv
import scanpy as sc
import velovi
import numpy as np
import pandas as pd

import os
import sys

import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler

import matplotlib.pyplot as plt
import seaborn as sns

import scvelo as scv

import anndata as ad

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *
from v4_functions import *

def get_adata(n_obs, mu, R, C, n_vars=1000, random_seed=317):
    cov = C.dot(C.T) * R
    alpha, beta, gamma = np.exp(np.random.multivariate_normal(mu, cov, size=n_vars).T)  # multivariate log-normal
    beta /= 3
    gamma /= 3
    # remove outliers
    idx = (alpha < np.percentile(alpha, 99)) & (beta < np.percentile(beta, 99)) & (gamma < np.percentile(gamma, 99))
    alpha = alpha[idx]
    beta = beta[idx]
    gamma = gamma[idx]
    n_vars = np.sum(idx)
    switches = np.random.uniform(.1, .5, size=n_vars)
    adata = scv.datasets.simulation(
        n_obs=n_obs,
        t_max=20,
        n_vars=n_vars,
        switches=switches,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        noise_level=.8,
        random_seed=random_seed
    )
    #scv.pp.neighbors(adata)
    #scv.tl.velocity(adata, mode='steady_state', use_raw=True)
    return adata

def run_countsplit_with_overdispersion(S,U,split_seed,overdisp_S,overdisp_U):
    print("########### Countsplitting")
    np.random.seed(split_seed)
    s1, s2  = countsplit(S,overdisps=overdisp_S)
    u1, u2  = countsplit(U,overdisps=overdisp_U)
    return [[s1,u1],[s2,u2]]

def create_adata(S_split,U_split,adata_total):
    adata_split = ad.AnnData(X=S_split.astype(np.float32))
    adata_split.layers["spliced"] = S_split
    adata_split.layers["unspliced"] = U_split
    adata_split.obs = pd.DataFrame(index=adata_total.obs.index)
    for obs_col in adata_total.obs.columns:
        adata_split.obs[obs_col] = adata_total.obs[obs_col].copy()
    adata_split.var = pd.DataFrame(index=adata_total.var.index)
    for var_col in adata_total.var.columns:
        adata_split.var[var_col] = adata_total.var[var_col].copy()
    return adata_split

def countsplit_and_create_adata(S,U,total,split_seed,overdisp_S,overdisp_U):
    print("########### Running the function for overdispersion estimation and countsplitting")
    split1,split2 = run_countsplit_with_overdispersion(S=S,U=U,split_seed=split_seed,overdisp_S=overdisp_S,overdisp_U=overdisp_U)
    print("########### Creating split adata objects")
    adata1 = create_adata(split1[0],split1[1],total)
    adata2 = create_adata(split2[0],split2[1],total)
    return adata1,adata2

mu = np.array([1, .2, .05])
R = np.array([[1., .2, .2],
              [.2, 1., .8],
              [.2, .8, 1.]])
C = np.array([.4, .4, .4])[:, None]

adata = get_adata(n_obs=1000, mu=mu, R=R, C=C, random_seed=317)
overdisp_S = estimate_overdisps(adata.layers['spliced'])
overdisp_U = estimate_overdisps(adata.layers['unspliced'])

adata_split1, adata_split2 = countsplit_and_create_adata(S=adata.layers['spliced'],
                                                         U=adata.layers['unspliced'],
                                                         total=adata,split_seed=317,
                                                         overdisp_S=overdisp_S,overdisp_U=overdisp_U)

adata.var['highly_variable'] = True
adata_split1.var['highly_variable'] = True
adata_split2.var['highly_variable'] = True

adata.obs['celltype'] = 'type1'
adata_split1.obs['celltype'] = 'type1'
adata_split2.obs['celltype'] = 'type1'

#sc.pp.neighbors(adata)
#sc.pp.neighbors(adata_split1)
#sc.pp.neighbors(adata_split2)

adata.write('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_total.h5ad')
adata_split1.write('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split1.h5ad')
adata_split2.write('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split2.h5ad')


def plot_cosine_similarity(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,text_x=None,text_y=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    print('median cos_sim='+str(np.median(cos_sim))+', mean cos_sim='+str(np.mean(cos_sim)))
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset+'_'+method
    # histogram
    print('Plot histogram')
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='gainsboro',color='powderblue') 
    max_frequency = np.max(counts)
    if text_x is None: text_x = np.quantile(cos_sim,[.05])[0]
    if text_y is None: text_y = max_frequency/5
    plt.axvline(np.mean(cos_sim), color='brown', linestyle='dashed', linewidth=1.5) ## add mean
    plt.axvline(np.median(cos_sim), color='peru', linestyle='dashed', linewidth=1.5) ## add median
    plt.text(text_x,text_y*2.5,'mean='+str(np.round(np.mean(cos_sim),4)), color='firebrick', fontsize=11)
    plt.text(text_x,text_y*3,'median='+str(np.round(np.median(cos_sim),4)), color='sienna', fontsize=11)
    plt.xlabel('cosine similarity')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, '+dataset+'+'+method)
    plt.savefig(fig_folder+dataset_method+'_cos_sim_hist.png')
    plt.clf()


## scv
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

total = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_total.h5ad')
split1 = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split1.h5ad')
split2 = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split2.h5ad')


def run_scv(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata,n_jobs=8) 
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata,n_jobs=8) 

run_scv(total)
run_scv(split1)
run_scv(split2)

import numpy as np
from scipy.stats import pearsonr,spearmanr

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import compute_cosine_similarity_union

fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2501-1/'
np.median(compute_cosine_similarity_union(split1,split2,'scv')[0])
plot_cosine_similarity(split1,split2,adata,dataset='sim',method='scv',
                       fig_folder='/root/capsule/results/',split_seed=317,
                       recompute=True,text_x=None,text_y=None)

plot_cosine_similarity(split1, split2, total, dataset='sim2501-1',method='scv',fig_folder=fig_folder)


arr_idx = np.where(~np.isnan(total.var['fit_gamma'] / total.var['fit_beta']))
pearsonr( (np.array(total.var['true_gamma'] / total.var['true_beta']))[arr_idx[0]], 
         (total.var['fit_gamma'] / total.var['fit_beta'])[arr_idx[0]] )

scv.tl.velocity_confidence(total)
np.median(total.obs['velocity_confidence'])

scv.tl.velocity_pseudotime(total, use_velocity_graph=False)
scv.tl.latent_time(total)
[spearmanr(total.obs.true_t, total.obs.velocity_pseudotime),
 spearmanr(total.obs.true_t, total.obs.latent_time)]

## utv
import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os

# the below script uses the environment: "utvClone"
# larry/erythroid configuration
velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True # default value
# (bool) linear regression $R^2$ on extreme quantile (default) or full data (adjusted)
# valid when self.VGENES = 'basic'
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1  # default initialization
# (float) threshold of R2 at later stage of the optimization proces
# to capture the dynamics of more genes beside initially selected velocity genes
# self.AGENES_R2 = 1 will switch to origianl mode with no amplification
os.environ["TF_USE_LEGACY_KERAS"]="1"

split1 = utv.run_model('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split1.h5ad', 
                       'celltype', config_file=velo_config)

split2 = utv.run_model('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split2.h5ad', 
                       'celltype', config_file=velo_config)

total = utv.run_model('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_total.h5ad', 
                       'celltype', config_file=velo_config)

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import compute_cosine_similarity_union


import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import compute_cosine_similarity_union
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2501-1/'
plot_cosine_similarity(split1, split2, total, dataset='sim2501-1',method='utv',fig_folder=fig_folder)
compute_cosine_similarity_union(split1,split2,'utv')

import numpy as np
from scipy.stats import pearsonr,spearmanr

arr_idx = np.where(~np.isnan(total.var['fit_gamma'] / total.var['fit_beta']))
pearsonr( (np.array(total.var['true_gamma'] / total.var['true_beta']))[arr_idx[0]], 
         (total.var['fit_gamma'] / total.var['fit_beta'])[arr_idx[0]] )

scv.tl.velocity_confidence(total)
np.median(total.obs['velocity_confidence'])

scv.tl.velocity_pseudotime(total, use_velocity_graph=False)
spearmanr(total.obs.true_t, total.obs.velocity_pseudotime)

## sct
split_seed = 317
method = 'sct'
sct_seed = 615

import sctour as sct
import scanpy as sc
import numpy as np
import random
import torch
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *

total = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_total.h5ad')
split1 = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split1.h5ad')
split2 = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split2.h5ad')

def sct_train_and_return_tnode(adata, sct_seed=615):
    import sctour as sct
    torch.manual_seed(sct_seed)
    random.seed(sct_seed)
    np.random.seed(sct_seed)
    adata.X = adata.X.astype(np.float32)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)
    #tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode = sct.train.Trainer(adata, loss_mode='mse', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
    return tnode

total = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_total.h5ad')

tnode_total = sct_train_and_return_tnode(total, sct_seed=sct_seed)
tnode_split1 = sct_train_and_return_tnode(split1, sct_seed=sct_seed)
tnode_split2 = sct_train_and_return_tnode(split2, sct_seed=sct_seed)

def compute_sct_avg_velocity(tnode,timesteps):
    v_shape = tnode.adata.shape
    v = np.zeros(v_shape)
    for t in timesteps:
        v += compute_sctour_velocity(tnode, timestep=t)
    return v/len(timesteps)

timesteps=[i/50 for i in range(1,11)]
total.layers['velocity'] = compute_sct_avg_velocity(tnode_total, timesteps)
split1.layers['velocity'] = compute_sct_avg_velocity(tnode_split1, timesteps) 
split2.layers['velocity'] = compute_sct_avg_velocity(tnode_split2, timesteps)

get_umap_sct(total)
get_umap_sct(split1)
get_umap_sct(split2)

## velovi
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import random
import datetime
from velovi import preprocess_data, VELOVI

def velovi_run_model(adata, method, split_seed):
    from velovi import preprocess_data, VELOVI
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if (not 'wop' in method): adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    return adata, vae
    
def add_velovi_outputs_to_adata(adata, vae, velovi_seed=615):
    torch.manual_seed(velovi_seed)
    random.seed(velovi_seed)
    np.random.seed(velovi_seed)
    latent_time = vae.get_latent_time(n_samples=25)
    velocities = vae.get_velocity(n_samples=25, velo_statistic="mean")
    t = latent_time
    scaling = 20 / t.max(0)
    adata.layers["velocity"] = velocities / scaling
    adata.layers["latent_time_velovi"] = latent_time
    adata.var["fit_alpha"] = vae.get_rates()["alpha"] / scaling
    adata.var["fit_beta"] = vae.get_rates()["beta"] / scaling
    adata.var["fit_gamma"] = vae.get_rates()["gamma"] / scaling
    adata.var["fit_t_"] = (
        torch.nn.functional.softplus(vae.module.switch_time_unconstr)
        .detach()
        .cpu()
        .numpy()
    ) * scaling
    adata.layers["fit_t"] = latent_time.values * np.array(scaling)[np.newaxis, :] # scaling[np.newaxis, :] 
    adata.var['fit_scaling'] = 1.0

total = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_total.h5ad')
split1 = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split1.h5ad')
split2 = sc.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2501-1_nobs1000_nvars1000_split2.h5ad')

split1, vae_split1 = velovi_run_model(split1, method='velovi', split_seed=317)
split2, vae_split2 = velovi_run_model(split2, method='velovi', split_seed=317)
total, vae_total = velovi_run_model(total, method='velovi', split_seed=317)

add_velovi_outputs_to_adata(split1, vae_split1, velovi_seed=615)
add_velovi_outputs_to_adata(split2, vae_split2, velovi_seed=615)
add_velovi_outputs_to_adata(total, vae_total, velovi_seed=615)

import numpy as np
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import compute_cosine_similarity_union
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2501-1/'
plot_cosine_similarity(split1, split2, total, dataset='sim2501-1',method='velovi',fig_folder=fig_folder)

import numpy as np
from scipy.stats import pearsonr,spearmanr

arr_idx = np.where(~np.isnan(total.var['fit_gamma'] / total.var['fit_beta']))
pearsonr( (np.array(total.var['true_gamma'] / total.var['true_beta']))[arr_idx[0]], 
         (total.var['fit_gamma'] / total.var['fit_beta'])[arr_idx[0]] )

scv.tl.velocity_confidence(total)
np.median(total.obs['velocity_confidence'])

scv.tl.velocity_graph(total)
scv.tl.velocity_pseudotime(total, use_velocity_graph=False)
spearmanr(total.obs.true_t, total.obs.velocity_pseudotime)
