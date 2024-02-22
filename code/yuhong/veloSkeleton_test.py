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
plt.savefig("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/test_1.png")

# create a graph layout same as the above UMAP positions
coords = [ [adata2000.obsm["X_umap"][:,0][nnx[i,0]], adata2000.obsm["X_umap"][:,1][nnx[i,0]]] for i in range(len(knots)) ]
layout_r = igraph.Layout(coords)

# Segment Skeleton
# leaving the graph as a whole in this case, but just run this for constructing the graph objects for the skeleton
import skeleton_methods.clustering as skelclus
knots = skeleton["centers"]
X_knotlabels = skeleton["cluster"]

#similarities
weights = skeleton["voron_weights"]
kcut = 1
hclustKnot = skelclus.cluster_weights(weights, X_knotlabels, kcut = kcut, method='average')

tol = 1e-10
#choose correct cut
cutheight = hclustKnot["hclust"][-kcut][2]
#cut some edges
weights[weights<(np.max(weights)-cutheight+tol)] = 0

#Euclidean distance between centers
# can be updated if doing cuts and distances between unconnected components are set to 0
kkdists = distance_matrix(knots,knots)
kkdists[weights == 0] = 0

skeleton.update({"kkdists": kkdists, "cutWeights": weights}) 

"""
import pickle
# save the computed skeleton
pickle.dump(skeleton, open("Msskeleton_r.pkl", "wb"))
# save the layout
pickle.dump(layout_r, open("Mslayout_r.pkl", "wb"))

# Just need to load with new code running process
import pickle
skeleton_r = pickle.load(open("Msskeleton_r.pkl","rb"))
layout_r = pickle.load(open("Mslayout_r.pkl","rb"))
skeleton_r.keys()
"""

# plotting the skeleton graph using the same layout us UMAP
g = igraph.Graph.Weighted_Adjacency(skeleton["kkdists"].tolist(),  mode='undirected')
fig, ax = plt.subplots()
# layout_r = g.layout("kk")
out = igraph.plot(g, layout=layout_r, target=ax)
plt.savefig("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/test_igraph_1.png")

# color the knots with the major cell type in each
from collections import defaultdict

knotClusbyCell = defaultdict(lambda: "Not Present")
for i in range(len(knots)):
    df = pd.Series(adata2000.obs['clusters'][X_nn[:,0] == i])
    knotClusbyCell[i] = df.value_counts()

majorCellClus = [knotClusbyCell[i].index[0] for i in range(len(knots))]
id_gen = igraph.UniqueIdGenerator()
color_indices = [id_gen.add(value) for value in majorCellClus]
palette = igraph.ClusterColoringPalette(len(id_gen))
colors = [palette[index] for index in color_indices]
g.vs["color"] = colors 
g.vs["label"] = [i for i in range(len(knots))]

fig, ax = plt.subplots(figsize=(12,8))
igraph.plot(g, layout=layout_r, target=ax, autocurve = True,edge_curved="0.1", 
            edge_width = 1,vertex_size = 0.5, vertex_label_size = 10)
plt.savefig("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/test_igraph_2.png")
