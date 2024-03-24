import scvelo as scv

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
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo_erythroid_v0.png")

scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)
scv.tl.score_genes_cell_cycle(adata)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)

from scipy.sparse import csr_matrix
adata.X = csr_matrix(adata.X)
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['Accession'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset

adata.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/adata_erythroid_v0.h5ad")


