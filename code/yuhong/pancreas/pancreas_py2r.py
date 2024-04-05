import scvelo as scv
import scanpy as sc
import scanpy.external as sce
from scipy.sparse import csr_matrix

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

adata = scv.datasets.pancreas()

spliced = adata.layers['spliced'].copy() # shape=(9815, 53801)
unspliced = adata.layers['unspliced'].copy()
gene_names = adata.var['highly_variable_genes'].copy()

scv.pp.filter_genes(adata, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) # only leave the top 2000 genes for this run for simplicity
scv.pp.log1p(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

print("Computing velocity")
# Using dynamical model (see paper) to calculate gene velocity
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata, n_jobs=8)

scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
scv.tl.score_genes_cell_cycle(adata)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)

adata.X = csr_matrix(adata.X)
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['highly_variable_genes'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset


adata.write_h5ad(filename = "~/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/adata_pancreas_preprocess.h5ad")

scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters',save="~/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/pancreas_umap_preprocess_clusters.png")
print("####### UMAPs plotted #######")


