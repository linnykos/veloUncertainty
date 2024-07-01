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

label = 'clusters'
label_color = 'clusters_colors'

## scatter plot without velocity arrows
scv.pl.scatter(total, color=label, cmap=label_color, save=figure_folder+"pancreas_preprocess_scatter.png")

## velocity confidence
scv.tl.velocity_confidence(total)
scv.pl.scatter(total, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], 
               save=figure_folder+"velo_conf/pancreas_preprocess_velo_confidence.png")

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
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    # ValueError: '_index' is a reserved name for dataframe columns.

add_velovi_outputs_to_adata(s1_317, vae_317s1)
add_velovi_outputs_to_adata(s2_317, vae_317s2)
# ValueError: Multi-dimensional indexing (e.g. `obj[:, None]`) is no longer supported. 
### Convert to a numpy array before indexing instead.

# write data
s1_317.write(filename=data_folder+"pancreas_seed317_split1.h5ad")
s2_317.write(filename=data_folder+"pancreas_seed317_split1.h5ad")
print("########################### adata splits saved ###########################")
vae_317s1.save(data_folder+'vae_seed317_split1.pt')
vae_317s2.save(data_folder+'vae_seed317_split2.pt')
print("########################### vae saved ###########################")

