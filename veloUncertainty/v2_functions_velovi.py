### conda activate veloviClone

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import random
from velovi import preprocess_data, VELOVI
import matplotlib.pyplot as plt
import seaborn as sns

#######################################
# add output
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

# compute umap
# for pancreas and larry
def compute_umap_pan(adata):
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # used to be 10
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)

def compute_umap_ery(adata):
    import bbknn
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("************ batch correction done ************")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # used to be 10
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)

#######################################
# plot velocity
def plot_velocity_velovi_pan(adata,adata_raw,dataset,method,fig_folder,fig_name,recompute=True): # the same as ery ver
    celltype_label = 'clusters'
    Ngenes = adata.layers['velocity'].shape[1]
    # umapCompute
    scv.pl.velocity_embedding_stream(adata, basis='umap', color=celltype_label,title="Velocity pan+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     recompute=recompute,save=fig_folder+"velocity/"+dataset+"_"+method+"_"+fig_name+"_umapCompute.png")
    # umapOriginal
    adata_plot = adata.copy()
    adata_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap', color=celltype_label,title="Velocity pan+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     recompute=recompute,save=fig_folder+"velocity/"+dataset+"_"+method+"_"+fig_name+"_umapOriginal.png")

def plot_velocity_velovi_ery(adata,adata_raw,dataset,method,fig_folder,fig_name,recompute=True):
    celltype_label = 'celltype'
    Ngenes = adata.layers['velocity'].shape[1]
    # umapCompute
    adata_plot = adata.copy()
    #scv.tl.velocity_graph(adata_plot)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap', color=celltype_label,title="Velocity ery+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     recompute=recompute,save=fig_folder+"velocity/"+dataset+"_"+method+"_"+fig_name+"_umapCompute.png")
    # umapOriginal
    adata_plot = adata.copy()
    adata_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    #scv.tl.velocity_graph(adata_plot)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap', color=celltype_label,title="Velocity ery+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     recompute=recompute,save=fig_folder+"velocity/"+dataset+"_"+method+"_"+fig_name+"_umapOriginal.png")

def plot_velocity_velovi(adata,adata_raw,dataset,method,fig_folder,fig_name,celltype_label=None,recompute=True):
    if dataset=='ery': celltype_label='celltype'
    elif 'pan' in dataset: celltype_label='clusters'
    elif dataset=='larry': celltype_label='state_info'
    Ngenes = adata.layers['velocity'].shape[1]
    # umapCompute
    adata_plot = adata.copy()
    scv.pl.velocity_embedding_stream(adata_plot,basis='umap',color=celltype_label,recompute=recompute,
                                     title='Velocity '+dataset+"+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     save=fig_folder+"velocity/"+dataset+"_"+method+"_"+fig_name+"_umapCompute.png")
    # umapOriginal
    adata_plot = adata.copy()
    adata_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    scv.pl.velocity_embedding_stream(adata_plot,basis='umap',color=celltype_label,recompute=recompute,
                                     title='Velocity '+dataset+"+velovi, "+fig_name+' (Ngenes='+str(Ngenes)+')',
                                     save=fig_folder+"velocity/"+dataset+"_"+method+"_"+fig_name+"_umapOriginal.png")


#######################################
# intrinsic uncertainty
def plot_uncertainty_velovi(adata, dataset, uncertainty_type, umap_type, fig_folder, fig_name, color_label=None,recompute=True,celltype_label=None):
    import os
    if dataset=='ery': celltype_label = 'celltype'
    elif 'pan' in dataset: celltype_label = 'clusters'
    elif dataset=='larry' and celltype_label==None: celltype_label = 'state_info'
    if 'C' in umap_type or 'c' in umap_type:
        umap_type = 'umapCompute'
    else: umap_type = 'umapOriginal'
    if 'in' in uncertainty_type or 'In' in uncertainty_type:
        uncertainty_type = 'intrinsic'
        color_label = 'directional_cosine_sim_variance'
        vmin = 'p0'
        vmax = 'p100'
    else: 
        uncertainty_type = 'extrinsic'
        color_label = 'directional_cosine_sim_variance_extrinisic'
        vmin = np.min(adata.obs[color_label])
        vmax = np.max(adata.obs[color_label])
        print(vmin,vmax)
    Ngenes = adata.layers['velocity'].shape[1]
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',recompute=recompute,
                                     title="Velocity "+dataset+'+velovi '+fig_name, frameon=False,size=100,alpha=0.5)
    sc.pl.umap(adata,color=color_label,cmap='coolwarm',ax=axs[1],legend_loc='none',vmin=vmin,vmax=vmax,
                title=uncertainty_type+' uncertainty, '+dataset+'+velovi '+fig_name+', Ngenes='+str(Ngenes), frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+'uncertainty/'+dataset+'_velovi_uncertainty_'+uncertainty_type+'_withRef_'+fig_name+'_'+umap_type+'.png')
    plt.clf()
    save_path = os.path.join(fig_folder, 'uncertainty/'+dataset+'_velovi_uncertainty_'+uncertainty_type+'_'+fig_name+'_'+umap_type+'.png')
    sc.pl.umap(adata, color=color_label, cmap="coolwarm", frameon=False, vmin=vmin,vmax=vmax,
               title=dataset+'+velovi, '+fig_name, show=False)
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()
    

def compute_intrinisic_uncertainty(adata_in,vae,adata_raw,dataset,fig_folder,fig_name,sample_seed=2216,n_samples=100,recompute=True,celltype_label=None):
    adata = adata_in.copy()
    import random
    torch.manual_seed(sample_seed)
    random.seed(sample_seed)
    np.random.seed(sample_seed)
    uncertainty_df, _ = vae.get_directional_uncertainty(n_samples=n_samples)
    uncertainty_df.head()
    for c in uncertainty_df.columns:
        adata.obs[c] = np.log10(uncertainty_df[c].values)
    # umapCompute
    plot_uncertainty_velovi(adata=adata,dataset=dataset,uncertainty_type='int',umap_type='umapCompute',
                            fig_folder=fig_folder,fig_name=fig_name,recompute=recompute,celltype_label=celltype_label)
    # umapOriginal
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap']
    plot_uncertainty_velovi(adata=adata,dataset=dataset,uncertainty_type='int',umap_type='umapOriginal',
                            fig_folder=fig_folder,fig_name=fig_name,recompute=recompute,celltype_label=celltype_label)

#######################################
# extrinsic uncertaity
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

def compute_extrinisic_uncertainty(adata_in,vae,adata_raw,dataset,fig_folder,fig_name,sample_seed=2216,n_samples=25,recompute=True,celltype_label=None):
    adata = adata_in.copy()
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
    plot_uncertainty_velovi(adata=adata,dataset=dataset,uncertainty_type='ext',umap_type='umapCompute',
                            fig_folder=fig_folder,fig_name=fig_name,recompute=recompute,celltype_label=celltype_label)
    # umapOriginal
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap']
    plot_uncertainty_velovi(adata=adata,dataset=dataset,uncertainty_type='ext',umap_type='umapOriginal',
                            fig_folder=fig_folder,fig_name=fig_name,recompute=recompute,celltype_label=celltype_label)

#######################################
# plot permutation score
def compute_permutation_score(adata,vae,dataset,fig_folder,fig_name):
    celltype_label = 'celltype'
    if 'pan' in dataset: celltype_label = 'clusters'
    perm_df, _ = vae.get_permutation_scores(labels_key=celltype_label)
    adata.var["permutation_score"] = perm_df.max(1).values
    plt.clf()
    sns.kdeplot(data=adata.var, x="permutation_score")
    plt.savefig(fig_folder+"uncertainty/"+dataset+'_velovi_permutation_score_'+fig_name+".png")
