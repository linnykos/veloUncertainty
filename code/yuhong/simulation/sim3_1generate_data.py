import scvelo as scv
import numpy as np
import scanpy as sc
import math
import matplotlib.pyplot as plt

def generate_negative_binomial_matrix(mean_matrix, overdispersion):
    # Ensure overdispersion is positive
    if overdispersion <= 0:
        raise ValueError("Overdispersion parameter must be positive.")
    
    # Calculate the parameters r and p for each element in the matrix
    r_matrix = np.full(mean_matrix.shape, overdispersion)
    p_matrix = r_matrix / (r_matrix + mean_matrix)
    
    # Generate the negative binomial random variables for each element in the matrix
    negative_binomial_matrix = np.random.negative_binomial(r_matrix, p_matrix)
    
    return negative_binomial_matrix

def generate_data(alpha, beta, gamma, max_thre, n_obs, n_vars,
                  noise_level, overdispersion, seed, t_max):
    # generate S with Gaussian noise
    adata = scv.datasets.simulation(n_obs=n_obs, n_vars=n_vars, t_max=t_max, alpha=alpha, beta=beta, gamma=gamma,
                                    noise_level=noise_level, random_seed=seed)
    S = adata.layers['spliced']
    U = adata.layers['unspliced']
    adata.layers['true_velocity'] = beta*U-gamma*S
    S_lam = math.e**(math.log(max_thre)*S/np.max(S)) - 1 
    U_lam = math.e**(math.log(max_thre)*U/np.max(U)) - 1
    np.random.seed(seed)
    S_poi = generate_negative_binomial_matrix(mean_matrix=S_lam,
                                             overdispersion=overdispersion)
    U_poi = generate_negative_binomial_matrix(mean_matrix=U_lam,
                                             overdispersion=overdispersion)
    adata.layers['spliced'] = S_poi.copy()
    adata.layers['unspliced'] = U_poi.copy()
    adata.X = S_poi.copy() # need to make sure X is the same as spliced matrix
    # List of quantiles you want to find
    quantiles = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    # Find the quantiles of the flattened array
    quantile_results = np.quantile(adata.X, quantiles)
    print(f"The quantiles {quantiles} of the matrix are: {quantile_results}")
    
    return adata

# good scv
n_obs = 1000
n_vars = 100
t_max = 25
alpha = 5
beta = .3
gamma = .5
noise_level = 2
seed = 520
max_thre = 100
overdispersion = 1000 # the larger means the more like a Poisson

data = generate_data(alpha=alpha, beta=beta, gamma=gamma, max_thre=max_thre, n_obs=n_obs, n_vars=n_vars,
                     noise_level=noise_level, overdispersion=overdispersion, seed=seed, t_max=t_max)
data.write("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3good.h5ad")

## True velocity
# https://github.com/theislab/scvelo/issues/1212 
scv.pp.log1p(data)
sc.pp.pca(data)
sc.pp.neighbors(data, n_pcs=30, n_neighbors=30)
scv.pp.moments(data, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(data,t_max=25)
scv.tl.velocity_graph(data,vkey="true_velocity")
scv.pl.velocity_embedding_grid(data, basis='pca', color="true_t", vkey="true_velocity", arrow_length=2, arrow_size=2, min_mass=10,
                               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3good_pca_vture_arrow.png")
scv.pl.velocity_embedding_stream(data, basis='pca', color="true_t", vkey="true_velocity",
                               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3good_pca_vture.png")         
del data

## run scv
data2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3good.h5ad")
scv.pp.log1p(data2)
sc.pp.pca(data2)
sc.pp.neighbors(data2, n_pcs=30, n_neighbors=30)
scv.pp.moments(data2, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(data2)
scv.tl.velocity(data2, mode='dynamical')
scv.tl.velocity_graph(data2)
scv.pl.velocity_embedding_grid(data2, basis='pca', color="true_t", arrow_length=2, arrow_size=2, min_mass=10,
                               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3good_pca_vcompute_arrow.png")
scv.pl.velocity_embedding_stream(data2, basis='pca', color="true_t",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3good_pca_vcompute.png")
scv.tl.velocity_confidence(data2)
scv.pl.scatter(data2, c='velocity_confidence', basis='pca',cmap='coolwarm', perc=[5, 95],
               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3good_pca_vcompute_scvconf.png")
# histogram
plt.hist(data2.obs['velocity_confidence'], bins=30, edgecolor='black')
mean_conf = np.mean(data2.obs['velocity_confidence'])
plt.axvline(mean_conf, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.85, 150, 'mean confidence = '+str(mean_conf), color='blue', fontsize=10)
## add labels and title
plt.xlabel('velocity_confidence')
plt.ylabel('Frequency')
plt.title('Histogram of velocity confidence, sim3 (good scv)')
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3good_hist_vconf.png')
plt.clf()
del data2

# bad scv
n_obs = 1000
n_vars = 20
t_max = 25
alpha = 4
beta = .2
gamma = .1
noise_level = 2
seed = 520
max_thre = 100
overdispersion = 1
data = generate_data(alpha=alpha, beta=beta, gamma=gamma, max_thre=max_thre, n_obs=n_obs, n_vars=n_vars,
                     noise_level=noise_level, overdispersion=overdispersion, seed=seed, t_max=t_max)
data.write("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3bad.h5ad")

## True velocity
scv.pp.log1p(data)
sc.pp.pca(data)
sc.pp.neighbors(data, n_pcs=30, n_neighbors=30)
scv.pp.moments(data, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(data,t_max=25)
scv.tl.velocity_graph(data,vkey="true_velocity")
scv.pl.velocity_embedding_grid(data, basis='pca', color="true_t", vkey="true_velocity", arrow_length=2, arrow_size=2, min_mass=10,
                               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_vture_arrow.png")
scv.pl.velocity_embedding_stream(data, basis='pca', color="true_t", vkey="true_velocity",
                               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_vture.png")         
del data

## run scv
data2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3bad.h5ad")
scv.pp.log1p(data2)
sc.pp.pca(data2)
sc.pp.neighbors(data2, n_pcs=30, n_neighbors=30)
scv.pp.moments(data2, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(data2)
scv.tl.velocity(data2, mode='dynamical')
scv.tl.velocity_graph(data2)
scv.pl.velocity_embedding_grid(data2, basis='pca', color="true_t", arrow_length=2, arrow_size=2, min_mass=10,
                               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_vcompute_arrow.png")
scv.pl.velocity_embedding_stream(data2, basis='pca', color="true_t",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_vcompute.png")
scv.tl.velocity_confidence(data2)
scv.pl.scatter(data2, c='velocity_confidence', basis='pca',cmap='coolwarm', perc=[5, 95],
               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_pca_vcompute_scvconf.png")
# histogram
plt.hist(data2.obs['velocity_confidence'], bins=30, edgecolor='black')
mean_conf = np.mean(data2.obs['velocity_confidence'])
plt.axvline(mean_conf, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.55, 150, 'mean confidence = '+str(mean_conf), color='blue', fontsize=10)
## add labels and title
plt.xlabel('velocity_confidence')
plt.ylabel('Frequency')
plt.title('Histogram of velocity confidence, sim3 (bad scv)')
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3bad_hist_vconf.png')
plt.clf()


