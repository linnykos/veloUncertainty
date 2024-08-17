### Did not run through all

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

import matplotlib.pyplot as plt
import seaborn as sns

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_velovi/"
s1_317_name = "pancreas_seed317_split1.h5ad"
s2_317_name = "pancreas_seed317_split2.h5ad"
total_name = "pancreas_velovi_preprocess.h5ad"

figure_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/velovi/"

total = sc.read(data_folder+total_name)

label = 'clusters'
label_color = 'clusters_colors'

## cosine similarity on preprocessed umap
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## velocity confidence
scv.tl.velocity_confidence(total)
scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], 
               save=figure_folder+"velo_conf/pancreas_preprocess_velo_confidence.png")

# 317 split1
s1_317 = sc.read(data_folder+"pancreas_seed317_split1.h5ad") # (3696, 351)
vae_317s1 = VELOVI.load(data_folder+'vae_seed317_split1.pt',s1_317)
## ValueError: Number of vars in `adata_target` not the same as source. Expected: 336 Received: 351
# 317 split2
s2_317 = sc.read(data_folder+"pancreas_seed317_split2.h5ad")
vae_317s2 = VELOVI.load(data_folder+'vae_seed317_split2.pt')

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

add_velovi_outputs_to_adata(s1_317, vae_317s1)
add_velovi_outputs_to_adata(s2_317, vae_317s2)
# ValueError: Multi-dimensional indexing (e.g. `obj[:, None]`) is no longer supported. 
### Convert to a numpy array before indexing instead.


def recover_umap_and_plot_velocity(data, fig_name):
    scv.tl.velocity_graph(data)
    scv.pl.velocity_embedding_stream(data, basis="pca",save=figure_folder+"velocity/pancreas_velovi_pca_"+fig_name+".png")
    # use the original umap
    data.obsm['X_umap'] = total.obsm['X_umap'].copy()
    scv.pl.velocity_embedding_stream(data, basis="umap",save=figure_folder+"velocity/pancreas_velovi_umapOriginal_"+fig_name+".png")
    del data.obsm['X_umap']
    # recover umap
    scv.pp.moments(data, n_pcs=30, n_neighbors=30)
    sc.tl.pca(data, svd_solver="arpack")
    sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
    sc.tl.umap(data)
    scv.pl.velocity_embedding_stream(data, basis="umap",save=figure_folder+"velocity/pancreas_velovi_umapRecover_"+fig_name+".png")
    scv.pl.scatter(data, color=label, cmap=label_color, basis="umap", save=figure_folder+"pancreas_velovi_scatter_umapRecover_"+fig_name+".png")

print("start recover_umap_and_plot_velocity")
recover_umap_and_plot_velocity(s1_317, "seed317_split1")
recover_umap_and_plot_velocity(s2_317, "seed317_split2")
print("finished recover_umap_and_plot_velocity")

# intrinsic uncertainty, figure output "uncertainty"
def compute_intrinsic_uncertatinty_and_plot(vae, adata, save_name):
    uncertainty_df, _ = vae.get_directional_uncertainty(n_samples=100)
    for c in uncertainty_df.columns:
        adata.obs[c] = np.log10(uncertainty_df[c].values)
    sc.pl.umap(adata, color="directional_cosine_sim_variance",cmap="Greys",vmin="p1",vmax="p99",save=save_name)

print("start compute_intrinsic_uncertatinty_and_plot")
compute_intrinsic_uncertatinty_and_plot(vae_317s1, s1_317, "pan_velovi_in_uncertainty_seed317_split1.png")
compute_intrinsic_uncertatinty_and_plot(vae_317s2, s2_317, "pan_velovi_in_uncertainty_seed317_split2.png")
print("finished compute_intrinsic_uncertatinty_and_plot")

# extrinsic uncertainty
def compute_extrinisic_uncertainty(adata, vae, n_samples=25) -> pd.DataFrame:
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
            t_mat = scv.utils.get_transition_matrix(
                adata, vkey=vkey, self_transitions=True, use_negative_cosines=True
            )
            extrapolated_cells = np.asarray(t_mat @ adata.layers["Ms"])
            extrapolated_cells_list.append(extrapolated_cells)
    extrapolated_cells = np.stack(extrapolated_cells_list)
    df = _compute_directional_statistics_tensor(extrapolated_cells, n_jobs=-1, n_cells=adata.n_obs)
    return df

def compute_extrinsic_uncertainty_and_plot(vae, adata, save_name):
    ext_uncertainty_df = compute_extrinisic_uncertainty(adata, vae)
    df = ext_uncertainty_df[0]
    for c in df.columns:
        adata.obs[c + "_extrinisic"] = np.log10(df[c].values)
    sc.pl.umap(adata, color="directional_cosine_sim_variance_extrinisic", vmin="p1", vmax="p99", save=save_name)

print("start compute_extrinsic_uncertainty_and_plot")
compute_extrinsic_uncertainty_and_plot(vae_317s1, s1_317, "pan_velovi_seed317_split1_ex_uncertainty.png")
compute_extrinsic_uncertainty_and_plot(vae_317s2, s2_317, "pan_velovi_seed317_split2_ex_uncertainty.png")
print("finished compute_extrinsic_uncertainty_and_plot")


# permutation score
def compute_permutation_score_and_plot(vae, adata, label, fig_name):
    perm_df, _ = vae.get_permutation_scores(labels_key=label)
    adata.var["permutation_score"] = perm_df.max(1).values
    plt.clf()
    sns.kdeplot(data=adata.var, x="permutation_score")
    plt.savefig(figure_folder+fig_name+".png")
    plt.clf()

print("start compute_permutation_score_and_plot")
compute_permutation_score_and_plot(vae_317s1, s1_317, "clusters", "pan_velovi_seed317_split1_permutation_score")
compute_permutation_score_and_plot(vae_317s2, s2_317, "clusters", "pan_velovi_seed317_split2_permutation_score")
print("finished compute_permutation_score_and_plot")

