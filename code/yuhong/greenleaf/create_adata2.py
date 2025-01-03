import scvelo as scv
import anndata as ad
import scanpy as sc
import pandas as pd
#import mygene
from scipy.sparse import csr_matrix, save_npz, load_npz

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/greenleaf/'
# Load the sparse matrix
count_matrix = load_npz(data_folder + "rna_sparse_counts.npz")
spliced_matrix = load_npz(data_folder + "rna_sparse_spliced_counts.npz")
unspliced_matrix = load_npz(data_folder + "rna_sparse_unspliced_counts.npz")

#    Assuming `rna_cell_metadata.txt` has rows of cells and columns for metadata,
#    and that one of these columns (likely the index) matches the cell names in `counts.columns`.
cell_metadata = pd.read_csv("/home/users/y2564li/kzlinlab/data/greanleaf_brain_multiome/rna_cell_metadata.txt", sep='\t', index_col=0)

# Load row and column names
gene_names = pd.read_csv(data_folder + "rna_gene_names.txt", header=None)
cell_names = pd.read_csv(data_folder + "rna_cell_names.txt", header=None)


spliced_gene_names = pd.read_csv(data_folder + "rna_spliced_gene_names.txt", header=None)
spliced_cell_names = pd.read_csv(data_folder + "rna_spliced_cell_names.txt", header=None)

unspliced_gene_names = pd.read_csv(data_folder + "rna_unspliced_gene_names.txt", header=None)
unspliced_cell_names = pd.read_csv(data_folder + "rna_unspliced_cell_names.txt", header=None)

# Extract gene and cell lists (as Series for convenience)
count_genes = gene_names[0]      # gene_names read from the count matrix file
count_cells = cell_names[0]

spliced_genes = spliced_gene_names[0]
spliced_cells = spliced_cell_names[0]

unspliced_genes = unspliced_gene_names[0]
unspliced_cells = unspliced_cell_names[0]

# Find common genes and cells across all three sets
common_genes = set(count_genes).intersection(spliced_genes, unspliced_genes)
common_cells = set(count_cells).intersection(spliced_cells, unspliced_cells)

# Convert back to lists (or Series) and preserve the original order of the count matrix
common_genes = count_genes[count_genes.isin(common_genes)]
common_cells = count_cells[count_cells.isin(common_cells)]
# We now have ordered lists of genes and cells that are common to all three matrices.

# Next, we need to determine the indices of these common genes and cells in each matrix.
# For each matrix, we know the original rows correspond to genes and columns to cells.

# Create a mapping from gene_name to original row index for each matrix:
count_gene_to_idx = pd.Series(range(count_genes.shape[0]), index=count_genes)
spliced_gene_to_idx = pd.Series(range(spliced_genes.shape[0]), index=spliced_genes)
unspliced_gene_to_idx = pd.Series(range(unspliced_genes.shape[0]), index=unspliced_genes)

# Similarly for cells:
count_cell_to_idx = pd.Series(range(count_cells.shape[0]), index=count_cells)
spliced_cell_to_idx = pd.Series(range(spliced_cells.shape[0]), index=spliced_cells)
unspliced_cell_to_idx = pd.Series(range(unspliced_cells.shape[0]), index=unspliced_cells)

# Get the integer indices for the common genes and cells for each matrix
count_gene_idx = count_gene_to_idx[common_genes].values
count_cell_idx = count_cell_to_idx[common_cells].values

spliced_gene_idx = spliced_gene_to_idx[common_genes].values
spliced_cell_idx = spliced_cell_to_idx[common_cells].values

unspliced_gene_idx = unspliced_gene_to_idx[common_genes].values
unspliced_cell_idx = unspliced_cell_to_idx[common_cells].values

# Now subset each matrix by common genes and cells
# Remember the shape of the matrices should be (n_genes, n_cells) 
# If they are not in that shape (e.g., cells in rows, genes in columns),
# adjust indexing accordingly.

count_matrix_sub = count_matrix[ count_gene_idx[:, None], count_cell_idx ]
spliced_matrix_sub = spliced_matrix[ spliced_gene_idx[:, None], spliced_cell_idx ]
unspliced_matrix_sub = unspliced_matrix[ unspliced_gene_idx[:, None], unspliced_cell_idx ]

# Now all three matrices correspond to the same sets of genes and cells in the same order.
# common_genes and common_cells hold the ordered names.

cell_metadata_sub = cell_metadata.loc[common_cells]

count_matrix_sub.shape


adata = ad.AnnData(
    X=count_matrix_sub.T,
    obs=cell_metadata_sub,
    var=pd.DataFrame(data={'indices': range(len(common_genes))}, index=common_genes)
)

# 6. Add layers for spliced and unspliced counts
adata.layers["spliced"] = spliced_matrix_sub.T
adata.layers["unspliced"] = unspliced_matrix_sub.T

# Read the cluster names file
cluster_df = pd.read_csv("/home/users/y2564li/kzlinlab/data/greanleaf_brain_multiome/scRNA_Colors.txt", sep='\t')
# Create a dictionary mapping Cluster.ID (e.g. c0) to Cluster.Name (e.g. nIPC/GluN1)
cluster_map = dict(zip(cluster_df['cluster'], cluster_df['potential name']))
# Now map these to a new column in adata.obs
adata.obs['cluster_name'] = adata.obs['seurat_clusters'].map(cluster_map)

umap_df = pd.read_csv("/home/users/y2564li/kzlinlab/data/greanleaf_brain_multiome/umapModel_rna_hft_uwot.csv", index_col=0)

umap_df.index = cell_names[0]

# Subset umap_df to these cells in the same order as they appear in adata
umap_sub = umap_df.loc[adata.obs_names.intersection(common_cells)]

# Convert to numpy array
adata.obsm["X_umap_greenleaf"] = umap_sub.values

#sc.pl.embedding(adata, basis="X_umap_greenleaf", color="cluster_name", 
#                save='/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_greenleaf/umap.png')

# Define the cell types to keep
cell_types_to_keep = [
    "cycling RG/IPC",
    "radial glia 1 (oRG)",
    "intermediate progenitor",
    "excitatory neuron 1",
    "excitatory neuron 2",
    "excitatory neuron 3",
    "excitatory neuron 6",
    "excitatory neuron 4",
    "excitatory neuron 5",
    "excitatory neuron 7"
]

# Subset the AnnData object
adata_subset = adata[adata.obs["cluster_name"].isin(cell_types_to_keep)].copy()
adata_subset.var.index = adata_subset.var.index.astype(str)
adata_subset.var.index.name = str(adata_subset.var.index.name)

#sc.pl.embedding(adata_subset, basis="X_umap_greenleaf", color="cluster_name",legend_loc="on data",
#                save='/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_greenleaf/umap_subset.png')

adata_subset.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/greenleaf_total_allgenes.h5ad')

adata.var.index = adata.var.index.astype(str)
adata.var.index.name = str(adata.var.index.name)
adata.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/greenleaf/greenleaf_full.h5ad')

