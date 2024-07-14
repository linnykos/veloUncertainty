import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

method = 'velovi'
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split1_v2.h5ad")
vae_split1 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_split1_v2.pt', split1)

split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split2_v2.h5ad")
vae_split2 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_split2_v2.pt', split2)

total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_total_v2.h5ad")
vae_total = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_total_v2.pt', total)

raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

#######################################
def add_velovi_outputs_to_adata(adata, vae):
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

print_message_with_time("############## Add velovi outputs to adata")

add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)

######################################################
## compute umap
def compute_umap_ery(adata):
    import bbknn
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("************ batch correction done ************")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)

print_message_with_time("############## Compute umap")

compute_umap_ery(split1)
compute_umap_ery(split2)
compute_umap_ery(total)

#######################################
## plot velocity
def plot_velocity_velovi_ery(adata,adata_raw,fig_name):
    celltype_label = 'celltype'
    Ngenes = adata.layers['velocity'].shape[1]
    # umapCompute
    adata_plot = adata.copy()
    #scv.tl.velocity_graph(adata_plot)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap', color=celltype_label,title="ery+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     save=fig_folder+"velocity/"+dataset_short+"_"+method+"_"+fig_name+"_umapCompute.png")
    # umapOriginal
    adata_plot = adata.copy()
    adata_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    #scv.tl.velocity_graph(adata_plot)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap', color=celltype_label,title="ery+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     save=fig_folder+"velocity/"+dataset_short+"_"+method+"_"+fig_name+"_umapOriginal.png")

print_message_with_time("############## Plot velocity")

plot_velocity_velovi_ery(adata=split1,adata_raw=raw,fig_name="split1")
plot_velocity_velovi_ery(adata=split2,adata_raw=raw,fig_name="split2")
plot_velocity_velovi_ery(adata=total,adata_raw=raw,fig_name="total")

#######################################
## plot cosine similarity
#compute_cosine_similarity(adata_split1=split1,adata_split2=split2,method='velovi')

#scv.tl.velocity_graph(split1)
#scv.tl.velocity_graph(split2)
#scv.tl.velocity_graph(total
print_message_with_time("############## Plot cosine similarity")

plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)

#######################################
## plot velo_conf
print_message_with_time("############## Plot velo_conf")

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)

#######################################
## plot ptime
# ???

#######################################
## intrinsic uncertainty
"""
save_path = os.path.join(save_dir, "Writeup6_eprs_mg_scvi_ePRS.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    adata2,
    color=["ePRS"],
    frameon=False,
    title="By ePRS status",
    size=5,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')
"""
def compute_intrinisic_uncertainty(adata_in,vae,adata_raw,dataset,fig_folder,fig_name,sample_seed=2216,n_samples=100):
    adata = adata_in.copy()
    celltype_label = 'celltype'
    if dataset == 'pan': celltype_label = 'clusters'
    Ngenes = adata.layers['velocity'].shape[1]
    import os
    import random
    torch.manual_seed(sample_seed)
    random.seed(sample_seed)
    np.random.seed(sample_seed)
    uncertainty_df, _ = vae.get_directional_uncertainty(n_samples=n_samples)
    uncertainty_df.head()
    for c in uncertainty_df.columns:
        adata.obs[c] = np.log10(uncertainty_df[c].values)
    # umapCompute
    save_path = os.path.join(fig_folder, 'uncertainty/'+dataset+'_velovi_uncertainty_intrinsic_'+fig_name+"_umapCompute.png")
    sc.pl.umap(adata, color="directional_cosine_sim_variance", cmap="coolwarm", perc=[1,100], 
               frameon=False,title=dataset+'+velovi, '+fig_name, show=False)
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     title="Velocity "+dataset+'+velovi '+fig_name, frameon=False,size=100,alpha=0.5)
    scv.pl.umap(adata,color='directional_cosine_sim_variance',cmap='coolwarm',perc=[1,100],ax=axs[1],legend_loc='none',
                title='Intrinsic uncertainty, '+dataset+'+velovi '+fig_name+', Ngenes='+str(Ngenes), frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+'uncertainty/'+dataset+'_velovi_uncertainty_intrinsic_withRef_'+fig_name+'_umapCompute.png')
    plt.clf()
    # umapOriginal
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap']
    save_path = os.path.join(fig_folder, 'uncertainty/'+dataset+'_velovi_uncertainty_intrinsic_'+fig_name+"_umapOriginal.png")
    sc.pl.umap(adata, color="directional_cosine_sim_variance", cmap="coolwarm", perc=[1,100], 
               frameon=False,title=dataset+'+velovi, '+fig_name, show=False)
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     title="Velocity "+dataset+'+velovi '+fig_name, frameon=False,size=100,alpha=0.5)
    scv.pl.umap(adata,color='directional_cosine_sim_variance',cmap='coolwarm',perc=[1,100],ax=axs[1],legend_loc='none',
                title='Intrinsic uncertainty, '+dataset+'+velovi '+fig_name+', Ngenes='+str(Ngenes), frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+'uncertainty/'+dataset+'_velovi_uncertainty_intrinsic_withRef_'+fig_name+'_umapOriginal.png')
    plt.clf()

print_message_with_time("############## Plot intrinsic uncertainty for split1")
compute_intrinisic_uncertainty(adata_in=split1,vae=vae_split1,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split1")

print_message_with_time("############## Plot intrinsic uncertainty for split2")
compute_intrinisic_uncertainty(adata_in=split2,vae=vae_split2,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split2")

print_message_with_time("############## Plot intrinsic uncertainty for total")
compute_intrinisic_uncertainty(adata_in=total,vae=vae_total,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="total")

#######################################
## extrinsic uncertainty
def compute_extrinisic_uncertainty_df(adata, vae, n_samples=25) -> pd.DataFrame:
    from velovi._model import _compute_directional_statistics_tensor
    from scvi.utils import track
    from contextlib import redirect_stdout
    import io
    extrapolated_cells_list = []
    for i in track(range(n_samples)):
        with io.StringIO() as buf, redirect_stdout(buf):
            vkey = "velocities_velovi_{i}".format(i=i)
            v = vae.get_velocity(n_samples=1, velo_statistic="mean")
            adata.layers[vkey] = v
            scv.tl.velocity_graph(adata, vkey=vkey, sqrt_transform=False, approx=True)
            t_mat = scv.utils.get_transition_matrix( adata, vkey=vkey, self_transitions=True, use_negative_cosines=True )
            extrapolated_cells = np.asarray(t_mat @ adata.layers["Ms"])
            extrapolated_cells_list.append(extrapolated_cells)
    extrapolated_cells = np.stack(extrapolated_cells_list)
    df = _compute_directional_statistics_tensor(extrapolated_cells, n_jobs=-1, n_cells=adata.n_obs)
    return df

def compute_extrinisic_uncertainty(adata_in,vae,adata_raw,dataset,fig_folder,fig_name,sample_seed=2216,n_samples=25):
    adata = adata_in.copy()
    celltype_label = 'celltype'
    if dataset == 'pan': celltype_label = 'clusters'
    Ngenes = adata.layers['velocity'].shape[1]
    import os
    import random
    torch.manual_seed(sample_seed)
    random.seed(sample_seed)
    np.random.seed(sample_seed)
    ext_uncertainty_df = compute_extrinisic_uncertainty_df(adata,vae,n_samples) 
    df = ext_uncertainty_df[0]
    for c in df.columns:
        adata.obs[c + "_extrinisic"] = np.log10(df[c].values)
    # umapCompute
    save_path = os.path.join(fig_folder, 'uncertainty/'+dataset+'_velovi_uncertainty_extrinsic_'+fig_name+"umapCompute.png")
    sc.pl.umap(adata, color="directional_cosine_sim_variance_extrinisic", perc=[1,100],
               frameon=False,title=dataset+'+velovi, '+fig_name, show=False)
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     title="Velocity "+dataset+'+velovi '+fig_name, frameon=False,size=100,alpha=0.5)
    scv.pl.umap(adata,color='directional_cosine_sim_variance_extrinisic',perc=[1,100],ax=axs[1],legend_loc='none',
                title='Extrinsic uncertainty, '+dataset+'+velovi '+fig_name+', Ngenes='+str(Ngenes), frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+'uncertainty/'+dataset+'_velovi_uncertainty_extrinsic_withRef_'+fig_name+'_umapCompute.png')
    plt.clf()
    # umapOriginal
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap']
    save_path = os.path.join(fig_folder, 'uncertainty/'+dataset+'_velovi_uncertainty_extrinsic_'+fig_name+"umapOriginal.png")
    sc.pl.umap(adata, color="directional_cosine_sim_variance_extrinisic", perc=[1,100],
               frameon=False,title=dataset+'+velovi, '+fig_name, show=False)
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     title="Velocity "+dataset+'+velovi '+fig_name, frameon=False,size=100,alpha=0.5)
    scv.pl.umap(adata,color='directional_cosine_sim_variance_extrinisic',perc=[1,100],ax=axs[1],legend_loc='none',
                title='Extrinsic uncertainty, '+dataset+'+velovi '+fig_name+', Ngenes='+str(Ngenes), frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+'uncertainty/'+dataset+'_velovi_uncertainty_extrinsic_withRef_'+fig_name+'_umapOriginal.png')
    plt.clf()

print_message_with_time("############## Plot extrinsic uncertainty for split1")
compute_extrinisic_uncertainty(adata_in=split1,vae=vae_split1,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split1")

print_message_with_time("############## Plot extrinsic uncertainty for split2")
compute_extrinisic_uncertainty(adata_in=split2,vae=vae_split2,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split2")

print_message_with_time("############## Plot extrinsic uncertainty for total")
compute_extrinisic_uncertainty(adata_in=total,vae=vae_total,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="total")

#######################################
## permutation score
def compute_permutation_score(adata,vae,dataset,fig_folder,fig_name):
    celltype_label = 'celltype'
    if dataset == 'pan': celltype_label = 'clusters'
    perm_df, _ = vae.get_permutation_scores(labels_key=celltype_label)
    adata.var["permutation_score"] = perm_df.max(1).values
    plt.clf()
    sns.kdeplot(data=adata.var, x="permutation_score")
    plt.savefig(fig_folder+"uncertainty/"+dataset+'_velovi_permutation_score_'+fig_name+".png")

print_message_with_time("############## Plot permutation score for split1")
compute_permutation_score(adata=split1,vae=vae_split1,dataset=dataset_short,fig_folder=fig_folder,fig_name="split1")

print_message_with_time("############## Plot permutation score for split2")
compute_permutation_score(adata=split2,vae=vae_split2,dataset=dataset_short,fig_folder=fig_folder,fig_name="split2")

print_message_with_time("############## Plot permutation score for split3")
compute_permutation_score(adata=total,vae=vae_total,dataset=dataset_short,fig_folder=fig_folder,fig_name="total")

print_message_with_time("############## All done")
