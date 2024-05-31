import scvelo as scv
import numpy as np

n_obs = 1000
n_vars = 500
t_max = 25
alpha = 5
beta = .3
gamma = .5
noise_level = 1
seed = 520
max_thre = 10

# generate S with Gaussian noise
data = scv.datasets.simulation(n_obs=n_obs,n_vars=n_vars,t_max=t_max,alpha=alpha,beta=beta,gamma=gamma,noise_level=1,random_seed=seed)
## obs: true_t
## var: true_t_, true_alpha, true_beta, true_gamma, true_scaling
# S = S*log(max_thre)/max(S)
# rpoisson(lambda=exp(S)) -> lambda=max_thre^{S/max(S)}
S = data.layers['spliced']
U = data.layers['unspliced']
S_lam = max_thre**(S/np.max(S))
U_lam = max_thre**(U/np.max(U))

np.random.seed(929)
S_poi = np.random.poisson(lam=[s for s in S_lam])
U_poi = np.random.poisson(lam=[u for u in U_lam])
data.layers['spliced_gaussian'] = S.copy()
data.layers['unspliced_gaussian'] = U.copy()
data.layers['spliced'] = S_poi.copy()
data.layers['unspliced'] = U_poi.copy()
data.X = S_poi.copy() # need to make sure X is the same as spliced matrix
data.layers['true_velocity'] = beta*U_poi-gamma*S_poi

data.write(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim1.h5ad")

### try the Gaussian version 
import scanpy as sc
n_obs = 1000
n_vars = 500
t_max = 25
alpha = 5
beta = .3
gamma = .5
noise_level = 1
seed = 520
max_thre = 10
data = scv.datasets.simulation(n_obs=n_obs,n_vars=n_vars,t_max=t_max,alpha=alpha,beta=beta,gamma=gamma,noise_level=1,random_seed=seed)
data.layers['true_velocity'] = beta * data.layers['spliced'] - gamma * data.layers['unspliced']
scv.pp.normalize_per_cell(data,enforce=True)
scv.pp.log1p(data)
scv.pp.moments(data, n_pcs=30, n_neighbors=30)
sc.tl.pca(data, svd_solver="arpack")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.umap(data,alpha=alpha,gamma=gamma)
scv.tl.recover_dynamics(data,t_max=25)
scv.tl.velocity_graph(data,vkey="true_velocity")
scv.pl.velocity_embedding_stream(data, basis='pca',vkey="true_velocity", save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1gaussian_pca_vtrue.png")
scv.pl.velocity_embedding_stream(data, basis='umap',vkey="true_velocity", save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1gaussian_umap_vtrue.png")


data = scv.datasets.simulation(n_obs=n_obs,n_vars=n_vars,t_max=t_max,alpha=alpha,beta=beta,gamma=gamma,noise_level=1,random_seed=seed)
scv.pp.normalize_per_cell(data,enforce=True)
scv.pp.log1p(data)
scv.pp.moments(data, n_pcs=30, n_neighbors=30)
sc.tl.pca(data, svd_solver="arpack")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.umap(data)
scv.tl.recover_dynamics(data)
scv.tl.velocity(data, mode="dynamical")
scv.tl.velocity_graph(data)
scv.pl.velocity_embedding_stream(data, basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1gaussian_pca_vcompute.png")
scv.pl.velocity_embedding_stream(data, basis='umap', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1gaussian_umap_vcompute.png")

###
# generate S without noise
data = scv.datasets.simulation(n_obs=n_obs,n_vars=n_vars,t_max=t_max,alpha=alpha,beta=beta,gamma=gamma,noise_level=0,random_seed=seed)
scv.pp.normalize_per_cell(data,enforce=True)
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
scv.pl.velocity_embedding_stream(data, basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1true_pca_vtrue.png")
scv.pl.velocity_embedding_stream(data, basis='umap', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1true_umap_vtrue.png")


data.write(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim1_true.h5ad")
