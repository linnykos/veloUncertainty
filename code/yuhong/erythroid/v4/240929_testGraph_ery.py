import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

dataset_long='erythroid'
dataset_short = 'ery'
method='velovi_woprep'
outputAdded = ''
outputAdded = '_outputAdded'
split1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
split2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')

A1=(split1.uns['velocity_graph']+split1.uns['velocity_graph_neg']).todense()
A2=(split2.uns['velocity_graph']+split2.uns['velocity_graph_neg']).todense()

A1 = split1.uns['velocity_graph'].todense()
A2 = split2.uns['velocity_graph'].todense()

A1 = split1.obsp['connectivities'].todense()
A2 = split2.obsp['connectivities'].todense()

### 
# Frobenius norm
np.linalg.norm(A1 - A2, 'fro')

### 
# Compute eigenvalues of both matrices
# Spectral Similarity (Eigenvalue Comparison)
eigenvalues_A1 = np.linalg.eigvals(A1)
eigenvalues_A2 = np.linalg.eigvals(A2)
np.linalg.norm(eigenvalues_A1 - eigenvalues_A2)

###
# cosine similarity
from sklearn.metrics.pairwise import cosine_similarity
# Flatten the adjacency matrices
A1_flat = A1.flatten().reshape(1, -1)
A2_flat = A2.flatten().reshape(1, -1)
# Calculate cosine similarity
cosine_similarity(A1_flat, A2_flat)

# Calculate weighted Jaccard similarity: similarity between the edge weights of the two graphs
# useful if your graphs are sparse and you're interested in comparing the magnitudes of corresponding edges
intersection = np.minimum(A1, A2).sum()
union = np.maximum(A1, A2).sum()
intersection / union

### 
# wasserstein distance: treat the adjacency matrices as distributions of edge weights and compute the Wasserstein distance
from scipy.stats import wasserstein_distance
# Flatten the adjacency matrices to 1D arrays
A1_flat = A1.flatten()
A2_flat = A2.flatten()
# Calculate Wasserstein distance
wasserstein_distance(A1_flat, A2_flat)

###
# graph edit distance
import networkx as nx
# Convert adjacency matrices to directed weighted NetworkX graphs
import scipy.sparse as sp
# Assuming A1 and A2 are sparse adjacency matrices (as numpy arrays)
A1_sparse = sp.csr_matrix(A1)  # Compressed Sparse Row matrix
A2_sparse = sp.csr_matrix(A2)
# Convert to NetworkX graphs using the sparse matrix representation
G1 = nx.from_scipy_sparse_array(A1_sparse, create_using=nx.DiGraph)
G2 = nx.from_scipy_sparse_array(A2_sparse, create_using=nx.DiGraph)
# Assign edge weights from the adjacency matrix
for u, v, d in G1.edges(data=True): d['weight'] = A1[u, v]
for u, v, d in G2.edges(data=True): d['weight'] = A2[u, v]
nx.graph_edit_distance(G1, G2)



################
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd


dataset_long='erythroid'
dataset_short = 'ery'
method='scv'
outputAdded = ''
outputAdded = '_outputAdded'
total = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4'+outputAdded+'.h5ad')
split1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
split2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')

# conda install -c conda-forge leidenalg
sc.tl.leiden(total,resolution=200)
# total.obs.leiden, uns['leiden']
split1.obs['leiden']= total.obs.leiden.copy()
split2.obs['leiden']= total.obs.leiden.copy()

clus = total.obs.leiden.cat.categories.values
print(clus)
vclus_split1 = np.zeros((len(clus), split1.shape[1]))
vclus_split2 = np.zeros((len(clus), split2.shape[1]))
for i in clus:
    cell_idx = np.where(total.obs.leiden == i)[0]
    v1 = split1.layers['velocity'][cell_idx,:].mean(0)
    v2 = split2.layers['velocity'][cell_idx,:].mean(0)
    vclus_split1[int(i)] = v1
    vclus_split2[int(i)] = v2

df1 = pd.DataFrame(vclus_split1, columns=split1.var.index)
df2 = pd.DataFrame(vclus_split2, columns=split2.var.index)
genes_intersect = np.intersect1d(split1.var.index[~np.isnan(vclus_split1[0])], split2.var.index[~np.isnan(vclus_split2[0])])
cos_sim = np.diag(cosine_similarity(df1[genes_intersect], df2[genes_intersect]))

np.quantile(cos_sim,[0.,.25,.5,.75,1.])
# velovi_woprep
# resolution=1000 (6300+): array([0.64446605, 0.93250802, 0.95929734, 0.97469183, 0.99362212])
# resolution=500 (4292): array([0.68233284, 0.94789891, 0.96802082, 0.97960063, 0.99419311])
# resolution=200: array([0.81443528, 0.96282195, 0.97587794, 0.98394937, 0.994257  ])

## velovi
# resolution=1000 (6391): array([-0.17444786,  0.50927629,  0.68815717,  0.78828742,  0.96099764])
# resolution=500 (4335): array([-0.15433898,  0.57497961,  0.74422688,  0.82052484,  0.93185173])
# resolution=200 (2513): array([-0.04604135,  0.65926113,  0.79830058,  0.85396911,  0.92891009])

## sct
# resolution=1000 (5910): array([-0.89324285,  0.66048444,  0.78787271,  0.87131178,  0.98475357])
# resolution=500 (2996): array([-0.89683729,  0.68205643,  0.79454771,  0.87686856,  0.97181476])
# resolution=200 (1405): array([-0.86399161,  0.69288195,  0.79360589,  0.87967353,  0.97118552])

## utv
# resolution=1000 (6295): array([-0.99001039,  0.99336349,  0.99672864,  0.99697969,  0.99711677])
# resolution=500 (3895): array([-0.97644671,  0.99372412,  0.99672241,  0.99697936,  0.9971146 ])
# resolution=200 (1957): array([-0.95345592,  0.99375566,  0.99669192,  0.99697877,  0.99709991])

## scv
# resolution=1000 (6372): array([-0.08351236,  0.4795517 ,  0.60135885,  0.68292487,  0.87969494])
# resolution=500 (4292): array([0.0290541 , 0.5522501 , 0.64678134, 0.71535299, 0.88541292])
# resolution=200 (2379): array([0.1707176 , 0.62606727, 0.69842645, 0.75419993, 0.89353241])
