import anndata as ad
import numpy as np
import scipy.sparse as sp
import mygene

def directed_laplacian_score(A, x):
    """
    Compute the normalized directed Laplacian score for a given vector x and adjacency matrix A.
    Assumes A is already a CSR sparse matrix!
    
    Args:
        A: scipy sparse matrix (n x n), non-negative, possibly asymmetric (already in CSR format).
        x: 1D numpy array (length n)
    
    Returns:
        Normalized directed Laplacian score (float)
    """
    # Do NOT call .tocsr() here anymore
    x_i_squared = A.multiply(x[:, None] ** 2).sum()
    x_j_squared = A.multiply(x[None, :] ** 2).sum()
    x_cross = A.multiply(np.outer(x, x)).sum()
    
    score = x_i_squared + x_j_squared - 2 * x_cross
    
    var_x = np.var(x)
    
    if var_x == 0:
        return np.nan
    normalized_score = score / var_x
    
    return normalized_score


# Load the AnnData object
adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/seed317/scv/adata_glf_scv_total_v4.h5ad")

# Step 1: Check all values are non-negative
velocity_graph = adata.uns['velocity_graph']

# Step 2: Normalize each row to sum to 1
# Normalize the velocity graph rows first (you already have this)
row_sums = np.array(velocity_graph.sum(axis=1)).flatten()
row_sums_safe = np.where(row_sums == 0, 1, row_sums)
row_scaling = sp.diags(1.0 / row_sums_safe)
velocity_graph_normalized = row_scaling.dot(velocity_graph)
velocity_graph_normalized = velocity_graph_normalized.tocsr()

# Now loop over ALL genes (all columns of adata)
n_genes = adata.n_vars  # This is 2000 genes in your case
scores = []

for gene_idx in range(10):
    print(f"Working on gene {gene_idx + 1} out of {n_genes}")
    
    # Get the gene expression vector x (for all cells, single gene)
    x = adata[:, gene_idx].X
    if sp.issparse(x):
        x = np.array(x.todense()).flatten()
    else:
        x = np.array(x).flatten()
    
    # Compute the normalized directed Laplacian score
    score = directed_laplacian_score(velocity_graph_normalized, x)
    scores.append(score)

# Convert list to numpy array
scores = np.array(scores)

print("Finished computing all scores!")


# Initialize
mg = mygene.MyGeneInfo()
# List of all 2000 gene IDs from your dataset
ensembl_ids = adata.var_names.tolist()  # Get all gene IDs from your AnnData object
# Query the mapping
out = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
# Build a mapping: Ensembl ID -> gene symbol
id_to_symbol = {entry['query']: entry.get('symbol', None) for entry in out}
# Print a few mappings to check
for k, v in list(id_to_symbol.items())[:10]:
    print(f"{k} → {v}")

# Count how many genes had no symbol
n_missing = sum(1 for symbol in id_to_symbol.values() if symbol is None)
# Total number of genes
n_total = len(id_to_symbol)
# Print summary
print(f"{n_missing} out of {n_total} Ensembl IDs did not have a gene symbol ({n_missing/n_total:.2%}).")


# Get gene names from adata
gene_names = adata.var_names.tolist()
# Map Ensembl IDs to gene symbols
gene_symbols = [id_to_symbol.get(gene, None) for gene in gene_names]
# Create a DataFrame
scores_df = pd.DataFrame({
    'gene': gene_names,
    'symbol': gene_symbols,
    'score': scores
})


# Save to CSV
output_path = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/Writeup15_greenleaf_gene_laplacian_scores.csv"
scores_df.to_csv(output_path, index=False)

print(f"Saved scores to {output_path}")