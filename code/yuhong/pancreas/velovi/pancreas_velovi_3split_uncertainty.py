import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

import matplotlib.pyplot as plt
import seaborn as sns

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_velovi/"
s1_317_name = "pancreas_seed317_split1_seurat.h5ad"
s2_317_name = "pancreas_seed317_split2_seurat.h5ad"
total_name = "pancreas_velovi_preprocess.h5ad"

figure_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/velovi/"

total = sc.read(data_folder+total_name)

## scatter plot without velocity arrows
scv.pl.scatter(total, color='clusters', cmap='clusters_colors', save=figure_folder+"pancreas_preprocess_scatter.png")

## velocity confidence
scv.tl.velocity_confidence(total)
scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], 
               save=figure_folder+"pancreas_preprocess_velo_confidence.png")

## cosine similarity on preprocessed umap
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 317 split1
s1_317 = sc.read(data_folder+s1_317_name)
scv.pp.filter_and_normalize(s1_317, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(s1_317, n_pcs=30, n_neighbors=30)
s1_317 = preprocess_data(s1_317) 
# train and apply model
VELOVI.setup_anndata(s1_317, spliced_layer="Ms", unspliced_layer="Mu")
vae_317s1 = VELOVI(s1_317)
vae_317s1.train()
# 317 split2
s2_317 = sc.read(data_folder+s2_317_name)
scv.pp.filter_and_normalize(s2_317, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(s2_317, n_pcs=30, n_neighbors=30)
s2_317 = preprocess_data(s2_317) 
# train and apply model
VELOVI.setup_anndata(s2_317, spliced_layer="Ms", unspliced_layer="Mu")
vae_317s2 = VELOVI(s2_317)
vae_317s2.train()

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

scv.tl.velocity_graph(s1_317)
scv.pl.velocity_embedding_stream(s1_317, basis='umap', save=figure_folder+"seed317_split1_umap.png")

scv.tl.velocity_graph(s2_317)
scv.pl.velocity_embedding_stream(s2_317, basis='umap', save=figure_folder+"seed317_split2_umap.png")

# write data
s1_317.write(filename=data_folder+"pancreas_seed317_split1.h5ad")
s2_317.write(filename=data_folder+"pancreas_seed317_split1.h5ad")
print("Wrote s1_317 and s2_317")
vae_317s1.save(data_folder+'vae_seed317_split1.pt')
vae_317s2.save(data_folder+'vae_seed317_split2.pt')
print("Saved vae")

# intrinsic uncertainty
def compute_intrinsic_uncertatinty_and_plot(vae, adata, save_name):
    uncertainty_df, _ = vae.get_directional_uncertainty(n_samples=100)
    for c in uncertainty_df.columns:
        adata.obs[c] = np.log10(uncertainty_df[c].values)
    sc.pl.umap(adata, color="directional_cosine_sim_variance",cmap="Greys",vmin="p1",vmax="p99",save=save_name)

compute_intrinsic_uncertatinty_and_plot(vae_317s1, s1_317, "pan_velovi_in_uncertainty_seed317_split1.png")
compute_intrinsic_uncertatinty_and_plot(vae_317s2, s2_317, "pan_velovi_in_uncertainty_seed317_split2.png")

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

compute_extrinsic_uncertainty_and_plot(vae_317s1, s1_317, "pan_velovi_seed317_split1_ex_uncertainty.png")
compute_extrinsic_uncertainty_and_plot(vae_317s2, s2_317, "pan_velovi_seed317_split2_ex_uncertainty.png")

# permutation score
def compute_permutation_score_and_plot(vae, adata, label, fig_name):
    perm_df, _ = vae.get_permutation_scores(labels_key=label)
    adata.var["permutation_score"] = perm_df.max(1).values
    plt.clf()
    sns.kdeplot(data=adata.var, x="permutation_score")
    plt.savefig(figure_folder+fig_name+".png")
    plt.clf()

compute_permutation_score_and_plot(vae_317s1, s1_317, "clusters", "pan_velovi_seed317_split1_permutation_score")
compute_permutation_score_and_plot(vae_317s2, s2_317, "clusters", "pan_velovi_seed317_split2_permutation_score")

# cosine similarity
