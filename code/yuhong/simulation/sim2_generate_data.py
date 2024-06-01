import scvelo as scv
import numpy as np
import scanpy as sc
import math

n_obs = 1000
n_vars = 500
t_max = 25
noise_level = 1
seed = 530
max_thre = 100
# default values
alpha = 5
beta = .5
gamma = .3

# generate S with Gaussian noise
data = scv.datasets.simulation(n_obs=n_obs,n_vars=n_vars,t_max=t_max,random_seed=seed)
## obs: true_t
## var: true_t_, true_alpha, true_beta, true_gamma, true_scaling
# S = max_thre*exp(S/max(S))
# rpoisson(lambda=exp(S)) -> lambda=max_thre^{S/max(S)}
S = data.layers['spliced']
U = data.layers['unspliced']
S_lam = math.e**(math.log(100)*S/np.max(S)) - 1 + 1e-5
U_lam = math.e**(math.log(100)*U/np.max(U)) - 1 + 1e-5

np.random.seed(2324)
S_poi = np.random.poisson(lam=[s for s in S_lam])
U_poi = np.random.poisson(lam=[u for u in U_lam])
data.layers['spliced_gaussian'] = S.copy()
data.layers['unspliced_gaussian'] = U.copy()
data.layers['spliced'] = S_poi.copy()
data.layers['unspliced'] = U_poi.copy()
data.X = S_poi.copy() # need to make sure X is the same as spliced matrix
data.layers['true_velocity'] = beta*U_poi - gamma*S_poi

data.write(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2.h5ad")

scv.pp.log1p(data)
scv.pp.moments(data, n_pcs=30, n_neighbors=30)
sc.tl.pca(data, svd_solver="arpack")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.umap(data,alpha=alpha,gamma=gamma)
scv.tl.recover_dynamics(data,t_max=25)
scv.tl.velocity_graph(data,vkey="true_velocity")
scv.pl.velocity_embedding_stream(data, basis='pca',color="true_t",vkey="true_velocity", save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2_pca_vtrue.png")
scv.pl.velocity_embedding_stream(data, basis='umap',color="true_t",vkey="true_velocity", save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2_umap_vtrue.png")


### try the Gaussian version 
data = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2.h5ad")
data.layers['spliced'] = data.layers['spliced_gaussian']
data.layers['unspliced'] = data.layers['unspliced_gaussian']
data.X = data.layers['spliced']
data.layers['true_velocity_gaussian'] = beta * data.layers['spliced'] - gamma * data.layers['unspliced']
scv.pp.log1p(data)
scv.pp.moments(data, n_pcs=30, n_neighbors=30)
sc.tl.pca(data, svd_solver="arpack")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.umap(data,alpha=alpha,gamma=gamma)
scv.tl.recover_dynamics(data,t_max=25)
scv.tl.velocity_graph(data,vkey="true_velocity_gaussian")
scv.pl.velocity_embedding_stream(data, basis='pca',color="true_t",vkey="true_velocity_gaussian", 
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2_pca_vtrueGaussian.png")
scv.pl.velocity_embedding_stream(data, basis='umap',color="true_t",vkey="true_velocity_gaussian", 
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2_umap_vtrueGaussian.png")

## recovered velocity
data = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2.h5ad")
scv.pp.log1p(data)
scv.pp.moments(data, n_pcs=30, n_neighbors=30)
sc.tl.pca(data, svd_solver="arpack")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.umap(data)
scv.tl.recover_dynamics(data)
scv.tl.velocity(data, mode="dynamical")
scv.tl.velocity_graph(data)
scv.pl.velocity_embedding_stream(data,color="true_t",basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2_pca_vcompute.png")
scv.pl.velocity_embedding_stream(data,color="true_t",basis='umap', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2_umap_vcompute.png")

###
# generate S without noise
data = scv.datasets.simulation(n_obs=n_obs,n_vars=n_vars,t_max=t_max,noise_level=0,random_seed=seed)
scv.pp.log1p(data)
scv.pp.moments(data, n_pcs=30, n_neighbors=30)
# WARNING: You seem to have 109 duplicate cells in your data. Consider removing these via pp.remove_duplicate_cells.
# WARNING: The neighbor graph has an unexpected format (e.g. computed outside scvelo) or is corrupted (e.g. due to subsetting). Consider recomputing with `pp.neighbors`.computing moments based on connectivities
data.layers['velocity'] = beta*data.layers['unspliced']-gamma*data.layers['spliced']
sc.tl.pca(data, svd_solver="arpack")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.umap(data)
scv.tl.recover_dynamics(data)
scv.tl.velocity_graph(data)
# ValueError: Your neighbor graph seems to be corrupted. Consider recomputing via pp.neighbors.
scv.pl.velocity_embedding_stream(data,color="true_t",basis='pca',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2noise0_pca_vtrue.png")
scv.pl.velocity_embedding_stream(data,color="true_t",basis='umap',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim2/sim2noise0_umap_vtrue.png")


data.write(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim2_noise0.h5ad")
