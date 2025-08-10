from scvelo.preprocessing.neighbors import get_neighs

n_neighs = (get_neighs(adata, "distances") > 0).sum(1)

# Convert matrix to a 1D array (if needed)
n_neighs_array = np.array(n_neighs).flatten()

# Get unique values and their counts
unique_values, counts = np.unique(n_neighs_array, return_counts=True)

# Combine results into a table-like format
value_counts = dict(zip(unique_values, counts))
print(value_counts)


########

sc.pp.neighbors(adata, n_neighbors=100, n_pcs=5) 
n_neighs = (get_neighs(adata, "distances") > 0).sum(1)

# Convert matrix to a 1D array (if needed)
n_neighs_array = np.array(n_neighs).flatten()

# Get unique values and their counts
unique_values, counts = np.unique(n_neighs_array, return_counts=True)

# Combine results into a table-like format
value_counts = dict(zip(unique_values, counts))
print(value_counts)

############

sc.tl.pca(adata, svd_solver="arpack")
scv.pp.moments(adata, n_pcs=5, n_neighbors=100)
sc.pp.neighbors(adata, n_neighbors=100, n_pcs=5)

#########

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from anndata import AnnData

def guaranteed_connected_neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_pcs: int = None,
    use_rep: str = None,
    metric: str = 'euclidean',
    key_added: str = 'neighbors',
):
    """
    Compute a neighbors graph ensuring that each cell has exactly n_neighbors connections,
    guaranteeing the graph is connected.
    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    n_neighbors : int
        Number of neighbors to connect each cell to.
    n_pcs : int
        Number of principal components to use. If None, use all components.
    use_rep : str
        Representation to use for computing distances.
    metric : str
        Distance metric to use (default: 'euclidean').
    key_added : str
        Key under `adata.uns` and `adata.obsp` to save the graph.
    """
    # Choose representation
    if use_rep is not None:
        X = adata.obsm[use_rep]
    elif 'X_pca' in adata.obsm and n_pcs is not None:
        X = adata.obsm['X_pca'][:, :n_pcs]
    else:
        X = adata.X
    # Use cKDTree for fast neighbor search
    tree = cKDTree(X)
    distances, indices = tree.query(X, k=n_neighbors + 1)  # Include self (k+1)  
    # Remove self-loops (first column is the point itself)
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    # Build distance and connectivity matrices
    n_obs = adata.n_obs
    row = np.repeat(np.arange(n_obs), n_neighbors)
    col = indices.flatten()
    data = distances.flatten()
    # Sparse matrix representation of the weighted graph
    distances_matrix = csr_matrix((data, (row, col)), shape=(n_obs, n_obs))
    connectivities_matrix = (distances_matrix > 0).astype(np.float32)
    # Save to adata
    adata.uns[key_added] = {
        'params': {
            'n_neighbors': n_neighbors,
            'metric': metric,
            'method': 'guaranteed_connected',
        }
    }
    adata.obsp['distances'] = distances_matrix
    adata.obsp['connectivities'] = connectivities_matrix
    print(f"Graph with guaranteed {n_neighbors} neighbors per cell computed and stored in `adata.obsp`.")

# Example usage
guaranteed_connected_neighbors(adata, n_neighbors=15, use_rep='X_pca')
sc.tl.umap(adata)
scv.tl.velocity_graph(adata)

###########