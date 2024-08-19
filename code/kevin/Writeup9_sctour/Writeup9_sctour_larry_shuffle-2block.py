import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import random
import anndata as ad
import datetime
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

sct_seed = 615

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

adata = sc.read_h5ad("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry.h5ad")

torch.manual_seed(sct_seed)
random.seed(sct_seed)
np.random.seed(sct_seed)

################

# Check if the layers exist before trying to delete them
layers_to_remove = ['spliced', 'unspliced', 'ambiguous', 'matrix']

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

adata.layers['counts'] = adata.X.copy()

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

# Step 1: Identify genes in cluster 0
genes_cluster = adata.var.index[adata.var['gene_cluster'] == 0]
# Step 2: Extract the data for these genes
data_cluster = adata[:, genes_cluster].X
# Ensure data is dense if it’s sparse
data_cluster_dense = data_cluster.toarray() if not isinstance(data_cluster, np.ndarray) else data_cluster
# Step 3: Shuffle the data across all cells (rows) for these genes
# We shuffle along axis=0 (rows, i.e., cells)
np.random.shuffle(data_cluster_dense)
# Step 4: Replace the original data with the shuffled data
adata[:, genes_cluster].X = data_cluster_dense

##

# Step 1: Identify genes in cluster 0
genes_cluster = adata.var.index[adata.var['gene_cluster'] == 1]
# Step 2: Extract the data for these genes
data_cluster = adata[:, genes_cluster].X
# Ensure data is dense if it’s sparse
data_cluster_dense = data_cluster.toarray() if not isinstance(data_cluster, np.ndarray) else data_cluster
# Step 3: Shuffle the data across all cells (rows) for these genes
# We shuffle along axis=0 (rows, i.e., cells)
np.random.shuffle(data_cluster_dense)
# Step 4: Replace the original data with the shuffled data
adata[:, genes_cluster].X = data_cluster_dense

# If you want to convert back to a sparse matrix, you can do so (optional)
adata.X = csr_matrix(adata.X)

################

tnode = sct.train.Trainer(adata, 
                          loss_mode='nb', 
                          alpha_recon_lec=0.5, 
                          alpha_recon_lode=0.5)
tnode.train()
adata.obs['ptime'] = tnode.get_time()
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, 
                                         alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, 
                                            adata.obsm['X_TNODE'])

adata.write("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/Writeup0_sctour_larry_shuffle-2block.h5ad")
tnode.save_model(save_dir="/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/", 
                 save_prefix="Writeup9_sctour_larry_tnode_shuffle-2block")
print_message_with_time("########### Total data wrote")

