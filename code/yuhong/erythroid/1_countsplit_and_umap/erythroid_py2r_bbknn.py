import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

adata = scv.datasets.gastrulation_erythroid()

# save the count matrices
spliced = adata.layers['spliced'].copy() # shape=(9815, 53801)
unspliced = adata.layers['unspliced'].copy()
gene_names = adata.var['Accession'].copy()

scv.pp.filter_genes(adata, min_shared_counts=20)
### Filtered out 47456 genes that are detected 20 counts (shared).
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
### Extracted 2000 highly variable genes.
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata, batch_key='sequencing.batch')
adata.X = adata.X.toarray()
bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata)
bbknn.bbknn(adata, batch_key='sequencing.batch')
print("Batch correction done!")

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_preprocess_bbknn_celltype.png")
scv.pl.velocity_embedding_stream(adata, basis='umap', color='sequencing.batch',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_preprocess_bbknn_seqbatch.png")
print("UMAPs plotted!")

scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)
scv.tl.score_genes_cell_cycle(adata)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)


adata.X = csr_matrix(adata.X)
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['Accession'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset

print("Writing h5ad file!")
adata.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/adata_erythroid_preprocess_bbknn.h5ad")


