import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI

import matplotlib.pyplot as plt
import seaborn as sns

# load and preprocess data
adata = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad")
spliced = adata.layers['spliced'].copy() # shape=(9815, 53801)
unspliced = adata.layers['unspliced'].copy()
gene_names = adata.var['highly_variable_genes'].copy()

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# adata = preprocess_data(adata) # 3696 × 1074

# train and apply model
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()

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

add_velovi_outputs_to_adata(adata, vae)
# ValueError: Multi-dimensional indexing (e.g. `obj[:, None]`) is no longer supported. 
### Convert to a numpy array before indexing instead.

scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/velovi/tut_preprocess_umap.png")

# write data
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['highly_variable_genes'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset
adata.write(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_velovi/pancreas_velovi_preprocess.h5ad")

# intrinsic uncertainty
uncertainty_df, _ = vae.get_directional_uncertainty(n_samples=100)
uncertainty_df.head()

for c in uncertainty_df.columns:
    adata.obs[c] = np.log10(uncertainty_df[c].values)
sc.pl.umap(
    adata, color="directional_cosine_sim_variance",
    cmap="Greys", vmin="p1", vmax="p99", save="/tut_uncertainty_int.png")
# /home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/velovi/

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
ext_uncertainty_df = compute_extrinisic_uncertainty(adata, vae) # [3696 rows x 5 columns]
df = ext_uncertainty_df[0]

for c in df.columns:
    adata.obs[c + "_extrinisic"] = np.log10(df[c].values)
sc.pl.umap(adata, color="directional_cosine_sim_variance_extrinisic", vmin="p1", vmax="p99",save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/velovi/tut_uncertainty_ext.png")

# permutation score
perm_df, _ = vae.get_permutation_scores(labels_key="clusters")
adata.var["permutation_score"] = perm_df.max(1).values

plt.clf()
sns.kdeplot(data=adata.var, x="permutation_score")
plt.savefig("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/velovi/tut_permutation_score.png")
plt.clf()


