import scvelo as scv
import importlib
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization
adata2000 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo-test.h5ad")

# Manual compute cosine similarity - to verify the methodology given by scvelo

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import skeleton_methods.skeleton as skel
import numpy as np
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import distance_matrix
from scipy.stats import gaussian_kde
import collections
from collections import defaultdict
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score
from scipy.spatial import distance_matrix
import igraph
from numpy.linalg import norm
import random
 
random.seed(123)

# specify the location features and velocity
# Need to run this for following analysis
X = adata2000.layers["Ms"]
V = adata2000.layers["velocity"]

subset = np.ones(adata2000.n_vars, bool)
subset &= np.array( adata2000.var["velocity_genes"].values, dtype = bool)
X = X[:, subset]
V = V[:, subset]
nans = np.isnan(np.sum(V, axis=0))
if np.any(nans):
    X = X[:, ~nans]
    V = V[:, ~nans]
V -= np.nanmean(V, axis =1)[:,None]

# Sekeleton constructed based on Ms features
import skeleton_methods.skeleton as skel
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import distance_matrix
from scipy.stats import gaussian_kde
import collections
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score
from scipy.spatial import distance_matrix
import igraph
from collections import defaultdict
from numpy.linalg import norm
import random
 
random.seed(123)

# define the knot construction function
def constructKnots(X, centers = None, labels = None, k = None, rep = 100):
    """
      construct knots using overfitting kMeans
  
      Parameters
      ----------
      X : np.array
          the data ndarray
      centers: np.array of the knots, can be provided
      labels: np.array of the knot label that each data point belongs to
      k: the number of knots
      rep: times of random initialization for kMeans
  
      Returns
      -------
      centers: ndarray of the constructed knots
      cluster: array of knot labels to which the data points belong to
      nknots: number of knots
      withinss: within cluster sum of squares for each Voronoi cell corresponding to the knots
    """
    n, d = X.shape
    #construct knots
    if centers is None and labels is None:
        # Overfitting k-means
        #setting the number of knots k
        if k is None:
            k = round(np.sqrt(n))
        km = KMeans(n_clusters = k, n_init = rep)
        km.fit(X)
        centers = km.cluster_centers_
        labels = km.labels_
    
    elif labels is None:#centers provided but not labels
        nbrs = NearestNeighbors(n_neighbors=1).fit(centers)
        labels = nbrs.kneighbors(X, return_distance=False)
        labels = np.array([sublist[0] for sublist in labels])
        k = len(centers)
          
    elif centers is None:#labels provided but not centers
        elements_count = collections.Counter(labels)
        k = len(elements_count.items())
        centers = np.array([[0.0 for col in range(d)] for row in range(k)])
        for key, value in elements_count.items():
            centers[key,] = np.mean(X[labels == key,], axis=0)
    else:
        k = len(centers)
    withinss = np.array([0.0]*k)
    for i in range(k):
        withinss[i] = np.sum(np.square(X[labels == i,:]-centers[i,:]))
      
    return {"centers":centers, "cluster":labels, "nknots":k, "withinss": withinss}

# construct knots
conKnots = constructKnots(X, k = 100)
skeleton = conKnots
edgeWeights = skel.skelWeights(X, conKnots, wedge = ['voronoi'])
skeleton.update(edgeWeights)
X_nn = skeleton["nn"]
knots = skeleton["centers"]
kkdists = distance_matrix(knots,knots)

# Compute the cluster sizes
knotSizes = np.zeros(100)
for i in range(100):
    knotSizes[i]= sum(X_nn[:,0] == i)

#filter out small clusters and reconstruct
knots = np.array(knots[knotSizes>=5,:])
conKnots = constructKnots(X, centers= knots)
skeleton = conKnots
edgeWeights = skel.skelWeights(X, conKnots, wedge = ['voronoi'])
skeleton.update(edgeWeights)

# recompute cluster sizes
X_nn = skeleton["nn"]
knots = skeleton["centers"]
kkdists = distance_matrix(knots,knots)
knotSizes = np.zeros(len(knots))
for i in range(len(knots)):
    knotSizes[i]= sum(X_nn[:,0] == i)

# compute the voronoi similarity between knots, which is counting the data points with the pair of knots as the two closest knots.
# These weights can be used for approximate Delaunnay triangulation
def VoronSimilarity(knots, X_nn, kkdists):
    m = knots.shape[0]
    voron_weights = np.array([[0.0 for col in range(m)] for row in range(m)])
    for i in range(m-1): #loop through knots pairs
        # center1 = knots[i]
        wi1 = np.where(X_nn[:,0] == i)[0]
        wi2 = np.where(X_nn[:,1] == i)[0]
        for j in range(i+1,m):
            # center2 = knots[j]
            wj1 = np.where(X_nn[:,0] == j)[0]
            wj2 = np.where(X_nn[:,1] == j)[0]
            #data point indices within 2nn neighborhood of knots i, j
            nn2ij = np.union1d(np.intersect1d(wi1, wj2), np.intersect1d(wi2, wj1))
        
            if len(nn2ij) < 1 :#not in Delaunnay Triangulation
                voron_weights[i,j] = 0.0
                continue
        
            # Euclidean distance between two centers
            d12 = kkdists[i,j]

            #compute the Voronoi density 
            voron_weights[i,j] = len(nn2ij)/d12
    return voron_weights

voron_weights = VoronSimilarity(knots, X_nn, kkdists)
skeleton.update( {"voron_weights": voron_weights + np.transpose((voron_weights))})

# To get major cell types within each voronoi cell
from collections import defaultdict
X_nn = skeleton["nn"]
knots = skeleton["centers"]

knotClusbyCell = defaultdict(lambda: "Not Present")
for i in range(len(knots)):
    df = pd.Series(adata2000.obs['clusters'][X_nn[:,0] == i])
    knotClusbyCell[i] = df.value_counts()

# knots to further split and refine. Picked the cells where the second major cell type 
# has counts more than half of the most common cell type and larger than 10 in counts.
# the refine indices are manually picked for now
# but can be automated with the rules 
refineId = [17,24,35,38,55,78,81,88,99]
for id in refineId:
    knots[id,:] = np.mean( X[np.logical_and( X_nn[:,0] == id , adata2000.obs['clusters'] == knotClusbyCell[id].index[0]),:], axis= 0)
    knots = np.vstack((knots, 
                       np.mean( X[np.logical_and( X_nn[:,0] == id , adata2000.obs['clusters'] != knotClusbyCell[id].index[0]),:], axis= 0)))

## plot
# Use closest observation X for each knot for UMAP plotting. Save time for UMAP calculation
nbrs = NearestNeighbors(n_neighbors=1).fit(X)
nnx = nbrs.kneighbors(skeleton["centers"], return_distance=False) 
# make knots larger sizes for plotting
sizes = np.array([1]*X.shape[0])
sizes[nnx[:,0]] = 10
colors = ['red' if i in nnx[:,0] else 'blue' for i in range(X.shape[0])]
fig, ax = plt.subplots()
ax.scatter(adata2000.obsm["X_umap"][:,0],adata2000.obsm["X_umap"][:,1], s=sizes ,c=colors)
for i in range(len(knots)):
    ax.annotate(str(i), (adata2000.obsm["X_umap"][:,0][nnx[i,0]], adata2000.obsm["X_umap"][:,1][nnx[i,0]]))
plt.savefig("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/test_refine_knots_0.png")


### 
# recalculate skeleton quantities based on refined knots
conKnots = constructKnots(X, centers= np.array(knots))
skeleton_r = conKnots
edgeWeights = skel.skelWeights(X, conKnots, wedge = ['voronoi'])
skeleton_r.update(edgeWeights)

# examine the refined skeleton
X_nn = skeleton_r["nn"]
knots = skeleton_r["centers"]
kkdists = distance_matrix(knots,knots)
knotSizes = np.zeros(len(knots))
for i in range(len(knots)):
    knotSizes[i]= sum(X_nn[:,0] == i)
min(knotSizes)
# calculate the voronoi weights with refined knots
voron_weights_r = VoronSimilarity(knots, X_nn, kkdists)
skeleton_r.update( {"voron_weights": voron_weights_r + np.transpose((voron_weights_r))})
###

# Use closest observation X for each knot for UMAP plotting. Save time for UMAP calculation
nbrs = NearestNeighbors(n_neighbors=1).fit(X)
nnx = nbrs.kneighbors(skeleton_r["centers"], return_distance=False) 
# make knots larger sizes for plotting
sizes = np.array([1]*X.shape[0])
sizes[nnx[:,0]] = 10
colors = ['red' if i in nnx[:,0] else 'blue' for i in range(X.shape[0])]
fig, ax = plt.subplots()
ax.scatter(adata2000.obsm["X_umap"][:,0],adata2000.obsm["X_umap"][:,1], s=sizes ,c=colors)
for i in range(len(knots)):
    ax.annotate(str(i), (adata2000.obsm["X_umap"][:,0][nnx[i,0]], adata2000.obsm["X_umap"][:,1][nnx[i,0]]))
plt.savefig("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/test_refine_knots_1.png")