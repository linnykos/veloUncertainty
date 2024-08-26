import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import anndata as ad
import datetime
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

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

# Step 4: Compute the principal components
pca = PCA()
pca_result = pca.fit_transform(correlation_matrix)
pca_result = pca_result[:,0:29]

# Step 5: Perform k-means clustering to group the 2000 genes into 2 clusters
kmeans = KMeans(n_clusters=2, random_state=sct_seed)
clusters = kmeans.fit_predict(pca_result)

# Step 6: Store the cluster labels in adata.var
adata.var['gene_cluster'] = clusters

################

# Shuffling

adata.X = adata.layers['counts'].copy()

# Ensure that the data matrices are dense if they are sparse
adata_X_dense = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
adata_spliced_dense = adata.layers['spliced'].toarray() if not isinstance(adata.layers['spliced'], np.ndarray) else adata.layers['spliced']
adata_unspliced_dense = adata.layers['unspliced'].toarray() if not isinstance(adata.layers['unspliced'], np.ndarray) else adata.layers['unspliced']

# Step 1: Identify genes in cluster 1 (keeping the genes in cluster 0 intact)
genes_cluster = adata.var.index[adata.var['gene_cluster'] == 0].tolist()
# Convert the gene names to their corresponding indices in adata.var
genes_cluster_indices = [adata.var.index.get_loc(gene) for gene in genes_cluster]
# Determine the number of elements to select (half of the total)
n_to_select = int(len(genes_cluster_indices) / 2)
# Randomly select half of the indices
selected_indices = np.random.choice(genes_cluster_indices, size=n_to_select, replace=False)
# Step 2: Shuffle the data across all cells (rows) for these genes
# We shuffle along axis=0 (rows, i.e., cells)
shuffle_order = np.random.permutation(adata_X_dense.shape[0])  # Generate a random shuffle order for the rows (cells)
# Step 3: Replace the original data with the shuffled data
for i in genes_cluster_indices:
    # Apply the shuffle order
    adata_X_dense[:, i] = adata_X_dense[shuffle_order, i]
    adata_spliced_dense[:, i] = adata_spliced_dense[shuffle_order, i]
    adata_unspliced_dense[:, i] = adata_unspliced_dense[shuffle_order, i]

# Step 1: Identify genes in cluster 1 (keeping the genes in cluster 0 intact)
genes_cluster = adata.var.index[adata.var['gene_cluster'] == 1].tolist()
# Convert the gene names to their corresponding indices in adata.var
genes_cluster_indices = [adata.var.index.get_loc(gene) for gene in genes_cluster]
# Determine the number of elements to select (half of the total)
n_to_select = int(len(genes_cluster_indices) / 2)
# Randomly select half of the indices
selected_indices = np.random.choice(genes_cluster_indices, size=n_to_select, replace=False)
# Step 2: Shuffle the data across all cells (rows) for these genes
# We shuffle along axis=0 (rows, i.e., cells)
shuffle_order = np.random.permutation(adata_X_dense.shape[0])  # Generate a random shuffle order for the rows (cells)
# Step 3: Replace the original data with the shuffled data
for i in genes_cluster_indices:
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

adata.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup10/Writeup10_larry_2block-ws.h5ad")
adata_split1.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup10/Writeup10_larry_2block-ws_split1.h5ad")
adata_split2.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup10/Writeup10_larry_2block-ws_split2.h5ad")
print_message_with_time("########### Finish splitting")

