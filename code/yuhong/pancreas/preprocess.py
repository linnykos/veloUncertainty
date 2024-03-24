import scvelo as scv
import importlib
import anndata

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

print("Downloading data")
# standard preprocessing steps from scVelo
adata2000 = scv.datasets.pancreas()

scv.pp.filter_genes(adata2000, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata2000)
scv.pp.filter_genes_dispersion(adata2000, n_top_genes=2000) # only leave the top 2000 genes for this run for simplicity
scv.pp.log1p(adata2000)
scv.pp.moments(adata2000, n_pcs=30, n_neighbors=30)

print("Computing velocity")
# Using dynamical model (see paper) to calculate gene velocity
scv.tl.recover_dynamics(adata2000, n_jobs=8)
scv.tl.velocity(adata2000, mode='dynamical')
scv.tl.velocity_graph(adata2000, n_jobs=8)

adata2000.write_h5ad(filename = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo-test.h5ad")

# to read in this file, run the following line:
# adata2000 = anndata.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo-test.h5ad')

print("Done ! :)")