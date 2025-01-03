import scvelo as scv
import anndata as ad
import scanpy as sc
import pandas as pd
#import mygene
from scipy.sparse import csr_matrix, save_npz, load_npz

in_data_folder = '/home/users/y2564li/kzlinlab/data/greanleaf_brain_multiome/'
out_data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/greenleaf/'

# 1. Load raw counts
#    Assuming `rna_counts.tsv` has a header row of cell names and the first column as gene names.
counts = pd.read_csv(in_data_folder + "rna_counts.tsv.gz", sep='\t', index_col=0)

# 2. Load cell metadata
#    Assuming `rna_cell_metadata.txt` has rows of cells and columns for metadata,
#    and that one of these columns (likely the index) matches the cell names in `counts.columns`.
cell_metadata = pd.read_csv(in_data_folder + "rna_cell_metadata.txt", sep='\t', index_col=0)


sparse_counts = csr_matrix(counts.values)

# Suppose sparse_counts is a csr_matrix
# Save it to a file
save_npz(out_data_folder + "rna_sparse_counts.npz", sparse_counts)

pd.DataFrame(counts.index).to_csv(out_data_folder + "rna_gene_names.txt", index=False, header=False)
pd.DataFrame(counts.columns).to_csv(out_data_folder + "rna_cell_names.txt", index=False, header=False)

# 3. Load spliced and unspliced counts
#    Similar structure to raw counts, with gene names as row index and cell names as column headers.
spliced = pd.read_csv(in_data_folder + "rna_spliced_counts.tsv.gz", sep='\t', index_col=0)

sparse_spliced_counts = csr_matrix(spliced.values)
save_npz(out_data_folder + "rna_sparse_spliced_counts.npz", sparse_spliced_counts)
pd.DataFrame(spliced.index).to_csv(out_data_folder + "rna_spliced_gene_names.txt", index=False, header=False)
pd.DataFrame(spliced.columns).to_csv(out_data_folder + "rna_spliced_cell_names.txt", index=False, header=False)

unspliced = pd.read_csv(in_data_folder + "rna_unspliced_counts.tsv.gz", sep='\t', index_col=0)
sparse_unspliced_counts = csr_matrix(unspliced.values)
save_npz(out_data_folder + "rna_sparse_unspliced_counts.npz", sparse_unspliced_counts)
pd.DataFrame(unspliced.index).to_csv(out_data_folder + "rna_unspliced_gene_names.txt", index=False, header=False)
pd.DataFrame(unspliced.columns).to_csv(out_data_folder + "rna_unspliced_cell_names.txt", index=False, header=False)

# 4. Ensure that cell order and gene order align across all matrices
#    It's crucial that `counts`, `spliced`, and `unspliced` have the same gene index order and cell columns.
#    Also ensure that `cell_metadata`'s index matches the cell column names.

