import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import random
import datetime
from velovi import preprocess_data, VELOVI
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import read_data_v4,read_raw_adata,print_message_with_time,get_celltype_label

def velovi_run_model(data_version,dataset_long,dataset_short,method,data_folder,split_seed):
    from velovi import preprocess_data, VELOVI
    adata = None
    print_message_with_time("#################### "+dataset_long+'_'+method+'_'+data_version+": Read data ")
    adata = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version,allgenes=True,outputAdded=False)
    gene_names = adata.var.index.copy()
    S_mat = adata.layers['spliced'].copy()
    U_mat = adata.layers['unspliced'].copy()
    print_message_with_time("#################### "+data_version+": Preprocess data ")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if (not 'wop' in method): adata = preprocess_data(adata)
    # train and apply model
    print_message_with_time("#################### "+data_version+": Train model ")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    # save vae
    savedata_folder = data_folder+'seed'+str(split_seed)+'/'+method+'/'
    print_message_with_time("#################### "+data_version+": Save vae ")
    vae.save(savedata_folder+'vae_'+dataset_short+'_'+method+'_'+data_version+'_v4.pt',overwrite=True)
    # save original counts
    positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
    positions = [positions_dict[gene] for gene in adata.var.index]
    adata.layers['spliced_original'] = S_mat[:,positions]
    adata.layers['unspliced_original'] = U_mat[:,positions]
    # write data
    print_message_with_time("#################### "+data_version+": Save adata (final version) ")
    adata.write(filename=savedata_folder+'adata_'+dataset_short+'_'+method+'_'+data_version+'_v4.h5ad')
    print_message_with_time("#################### "+data_version+": All done for "+dataset_short+'_'+method+'_'+data_version)


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
def compute_umap(adata, dataset):
    if 'ery' in dataset: compute_umap_ery(adata)
    else:
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40) # used to be 10
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
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40) # used to be 10
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)

#######################################
import matplotlib.pyplot as plt
import seaborn as sns

#######################################
# intrinsic uncertainty
def plot_uncertainty_velovi(adata, dataset, uncertainty_type, umap_type, fig_folder, data_version, color_label=None,recompute=True,celltype_label=None):
    import os
    celltype_label = get_celltype_label(dataset)
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
                                     title="Velocity "+dataset+'+velovi '+data_version, frameon=False,size=100,alpha=0.5)
    sc.pl.umap(adata,color=color_label,cmap='coolwarm',ax=axs[1],legend_loc='none',vmin=vmin,vmax=vmax,
                title=uncertainty_type+' uncertainty, '+dataset+'+velovi '+data_version+', Ngenes='+str(Ngenes), frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+'uncertainty/'+dataset+'_velovi_uncertainty_'+uncertainty_type+'_withRef_'+data_version+'_'+umap_type+'.png')
    plt.clf()
    save_path = os.path.join(fig_folder, 'uncertainty/'+dataset+'_velovi_uncertainty_'+uncertainty_type+'_'+data_version+'_'+umap_type+'.png')
    sc.pl.umap(adata, color=color_label, cmap="coolwarm", frameon=False, vmin=vmin,vmax=vmax,
               title=dataset+'+velovi, '+data_version, show=False)
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()
    

def compute_intrinisic_uncertainty(adata_in,vae,dataset,fig_folder,data_version,sample_seed=2216,n_samples=100,recompute=True,celltype_label=None):
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
                            fig_folder=fig_folder,data_version=data_version,recompute=recompute,celltype_label=celltype_label)
    # umapOriginal
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal']
    plot_uncertainty_velovi(adata=adata,dataset=dataset,uncertainty_type='int',umap_type='umapOriginal',
                            fig_folder=fig_folder,data_version=data_version,recompute=recompute,celltype_label=celltype_label)

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

def compute_extrinisic_uncertainty(adata_in,vae,dataset,fig_folder,data_version,sample_seed=2216,n_samples=25,recompute=True,celltype_label=None):
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
                            fig_folder=fig_folder,data_version=data_version,recompute=recompute,celltype_label=celltype_label)
    # umapOriginal
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal']
    plot_uncertainty_velovi(adata=adata,dataset=dataset,uncertainty_type='ext',umap_type='umapOriginal',
                            fig_folder=fig_folder,data_version=data_version,recompute=recompute,celltype_label=celltype_label)

#######################################
# plot permutation score
def compute_permutation_score(adata,vae,dataset,fig_folder,data_version):
    celltype_label = 'celltype'
    if 'pan' in dataset: celltype_label = 'clusters'
    perm_df, _ = vae.get_permutation_scores(labels_key=celltype_label)
    adata.var["permutation_score"] = perm_df.max(1).values
    plt.clf()
    sns.kdeplot(data=adata.var, x="permutation_score")
    plt.savefig(fig_folder+"uncertainty/"+dataset+'_velovi_permutation_score_'+data_version+".png")