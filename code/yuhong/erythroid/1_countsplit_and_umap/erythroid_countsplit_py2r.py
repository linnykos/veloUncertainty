import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

adata = scv.datasets.gastrulation_erythroid()
# adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/data/Gastrulation/erythroid_lineage.h5ad")

spliced = adata.layers['spliced'].copy()
unspliced = adata.layers['unspliced'].copy()
gene_names = adata.var['Accession'].copy()

scv.pp.filter_genes(adata, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) # only leave the top 2000 genes for this run for simplicity
scv.pp.log1p(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

adata.X = csr_matrix(adata.X)
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['Accession'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset

adata.write_h5ad(filename = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_countsplit_py2r.h5ad")