import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
import datetime
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from scipy.sparse.linalg import eigsh

import sys
sys.path.append('/home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from countsplit import *

sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

adata = sc.read_h5ad("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry.h5ad")

random.seed(sct_seed)
np.random.seed(sct_seed)

################

# Check if the layers exist before trying to delete them
layers_to_remove = ['ambiguous', 'matrix']

for layer in layers_to_remove:
    if layer in adata.layers:
        del adata.layers[layer]

# Verify the layers have been removed
print(adata.layers)

################

adata.X = adata.X.astype(np.float32)
sc.pp.calculate_qc_metrics(adata, 
                           percent_top=None, 
                           log1p=False, 
                           inplace=True)
sc.pp.highly_variable_genes(adata, 
                            flavor='seurat_v3', 
                            n_top_genes=2000, 
                            subset=True)
##################

# Step 1: Log-transform the data
sc.pp.log1p(adata)  # This will log-transform adata.X in place

# Step 2: Compute the correlation matrix among all 2000 genes
# Convert to dense if adata.X is still in sparse format
data_dense = adata.X.toarray() 

# Check for NaN or infinite values and replace them with zeros or appropriate values
data_dense = np.nan_to_num(data_dense)

# Step 3: Compute the correlation matrix among all 2000 genes
correlation_matrix = np.corrcoef(data_dense.T)

# Step 4: Compute the eigendecomposition
# Compute the top 30 eigenvalues and corresponding eigenvectors
# Note: If there are fewer than 30 dimensions, adjust accordingly
k = min(30, correlation_matrix.shape[0])
eigenvalues, eigenvectors = eigsh(correlation_matrix, k=k, which='LM')

# # Plot the eigenvalues as a bar plot
# plt.figure(figsize=(10, 6))
# plt.bar(range(1, len(eigenvalues) + 1), eigenvalues, alpha=0.7, align='center')
# plt.xlabel('Eigenvalue Index')
# plt.ylabel('Eigenvalue')
# plt.title('Top 30 Eigenvalues')
# plt.grid(True)
# plt.show()

# Step 5: Cluster genes based on the eigenvectors corresponding to the top 5 eigenvalues
# (5 chosen here based on the barplots)
# Eigenvectors' columns correspond to the variables
top_eigenvectors = eigenvectors[:, :5]
# Perform K-means clustering to cluster the variables into 5 clusters
kmeans = KMeans(n_clusters=5, random_state=sct_seed).fit(top_eigenvectors)
# Get the cluster labels for each variable
cluster_labels = kmeans.labels_
# Step 6: Store the cluster labels in adata.var
adata.var['gene_cluster'] = cluster_labels
adata.var['gene_cluster'].value_counts()

################

# Shuffling

adata.X = adata.layers['counts'].copy()

# Ensure that the data matrices are dense if they are sparse
adata_X_dense = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
adata_spliced_dense = adata.layers['spliced'].toarray() if not isinstance(adata.layers['spliced'], np.ndarray) else adata.layers['spliced']
adata_unspliced_dense = adata.layers['unspliced'].toarray() if not isinstance(adata.layers['unspliced'], np.ndarray) else adata.layers['unspliced']

# Loop over all clusters
for cluster_label in adata.var['gene_cluster'].unique():
    print(f"Working on cluster {cluster_label}")
    # Step 1: Identify genes in the current cluster
    genes_cluster = adata.var.index[adata.var['gene_cluster'] == cluster_label].tolist()
    # Convert the gene names to their corresponding indices in adata.var
    genes_cluster_indices = [adata.var.index.get_loc(gene) for gene in genes_cluster]
    
    # Randomly shuffle the list of gene indices
    np.random.shuffle(genes_cluster_indices)
    
    # Step 2: Shuffle the data for both randomly selected halves of the genes
    for half in [0, 1]:
        # Determine the range for the current half
        start_idx = half * int(len(genes_cluster_indices) / 2)
        end_idx = (half + 1) * int(len(genes_cluster_indices) / 2)
        selected_indices = genes_cluster_indices[start_idx:end_idx]
         
        # Shuffle the data across all cells (rows) for these genes
        shuffle_order = np.random.permutation(adata_X_dense.shape[0])  # Generate a random shuffle order for the rows (cells)
        
        # Step 3: Replace the original data with the shuffled data
        for i in selected_indices:
            # Apply the shuffle order
            adata_X_dense[:, i] = adata_X_dense[shuffle_order, i]
            adata_spliced_dense[:, i] = adata_spliced_dense[shuffle_order, i]
            adata_unspliced_dense[:, i] = adata_unspliced_dense[shuffle_order, i]


# Replace the original data with the shuffled data
adata.X = csr_matrix(adata_X_dense)
adata.layers['spliced'] = csr_matrix(adata_spliced_dense)
adata.layers['unspliced'] = csr_matrix(adata_unspliced_dense)

################

# now split
# code from: https://github.com/linnykos/veloUncertainty/blob/main/code/yuhong/larry/allgenes_countsplit.py

S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

overdisps_S = estimate_overdisps(S_mat)
overdisps_U = estimate_overdisps(U_mat)

# Define the percentiles you want to compute (0%, 10%, 20%, ..., 100%)
percentiles = np.arange(0, 101, 10)
overdisps_S = np.nan_to_num(overdisps_S)
overdisps_U = np.nan_to_num(overdisps_U)
# Compute the quantiles at the specified percentiles
np.percentile(overdisps_S, percentiles)
np.percentile(overdisps_U, percentiles)

# Truncate values in the vector so that no value exceeds 1
# This is a janky fix
overdisps_S = np.minimum(overdisps_S, 1)
overdisps_U = np.minimum(overdisps_U, 1)

np.random.seed(sct_seed)
S_split1, S_split2  = countsplit(S_mat,overdisps=overdisps_S)
np.random.seed(sct_seed)
U_split1, U_split2  = countsplit(U_mat,overdisps=overdisps_U)

adata_split1 = ad.AnnData(X=S_split1.astype(np.float32))
adata_split1.layers["spliced"] = S_split1
adata_split1.layers["unspliced"] = U_split1
adata_split1.obs = pd.DataFrame(index=adata.obs.index)
for obs_col in adata.obs.columns:
    adata_split1.obs[obs_col] = adata.obs[obs_col].copy()

adata_split1.var = pd.DataFrame(index=adata.var.index)
for var_col in adata.var.columns:
    adata_split1.var[var_col] = adata.var[var_col].copy()

###

adata_split2 = ad.AnnData(X=S_split2.astype(np.float32))
adata_split2.layers["spliced"] = S_split2
adata_split2.layers["unspliced"] = U_split2
adata_split2.obs = pd.DataFrame(index=adata.obs.index)
for obs_col in adata.obs.columns:
    adata_split2.obs[obs_col] = adata.obs[obs_col].copy()

adata_split2.var = pd.DataFrame(index=adata.var.index)
for var_col in adata.var.columns:
    adata_split2.var[var_col] = adata.var[var_col].copy()

adata.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/Writeup11_larry_full-block.h5ad")
adata_split1.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/Writeup11_larry_full-block_split1.h5ad")
adata_split2.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/Writeup11_larry_full-block_split2.h5ad")
print_message_with_time("########### Finish splitting")

