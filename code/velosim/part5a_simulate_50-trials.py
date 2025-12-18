# scvelo virtual environment

import scvelo as scv
import numpy as np
import scanpy as sc
import scipy.io
import scipy.sparse
from pathlib import Path

def generate_data(alpha, beta, gamma, n_obs, n_vars, noise_level, seed, t_max):
    adata = scv.datasets.simulation(
        n_obs=n_obs,
        n_vars=n_vars,
        t_max=t_max,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        noise_level=noise_level,
        random_seed=seed,
    )
    
    scale_factor = 100.0
    
    # Ensure dense arrays before np.maximum (in case layers are sparse)
    S = np.asarray(adata.layers["spliced"])
    U = np.asarray(adata.layers["unspliced"])
    
    S_scaled = np.maximum(S * scale_factor, 0)
    U_scaled = np.maximum(U * scale_factor, 0)
    
    adata.layers["true_velocity"] = (beta * U_scaled) - (gamma * S_scaled)
    
    np.random.seed(seed)
    S_poisson = np.random.poisson(S_scaled).astype(np.float32)
    U_poisson = np.random.poisson(U_scaled).astype(np.float32)
    
    adata.layers["spliced"] = S_poisson
    adata.layers["unspliced"] = U_poisson
    adata.layers["spliced_raw"] = S_poisson.copy()
    adata.layers["unspliced_raw"] = U_poisson.copy()
    
    # Scanpy uses .X for PCA by default
    adata.X = S_poisson
    
    # PCA
    sc.pp.pca(adata)
    n_pcs_available = adata.obsm["X_pca"].shape[1]
    n_pcs_use = min(30, n_pcs_available)
    
    # Neighbors
    sc.pp.neighbors(adata, n_pcs=n_pcs_use, n_neighbors=30)
    
    # Moments
    scv.pp.moments(adata)
    
    quantiles = [0, 0.1, 0.5, 0.9, 1]
    quantile_results = np.quantile(adata.layers["spliced_raw"], quantiles)
    print(f"The quantiles {quantiles} of the new Poisson matrix are: {quantile_results}")
    
    return adata


# --- Usage ---
n_obs = 1000
n_vars = 50
t_max = 25
alpha = 5
beta = 0.3
gamma = 0.5
noise_level = 0.1

out_dir = Path("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/simulation")
out_dir.mkdir(parents=True, exist_ok=True)

for seed in range(50):
    print(f"Seed: {seed}")
    
    data = generate_data(
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        n_obs=n_obs,
        n_vars=n_vars,
        noise_level=noise_level,
        seed=seed,
        t_max=t_max,
    )
    
    spliced_sparse = scipy.sparse.csc_matrix(data.layers["spliced_raw"])
    unspliced_sparse = scipy.sparse.csc_matrix(data.layers["unspliced_raw"])
    
    scipy.io.mmwrite(out_dir / f"spliced_seed{seed}.mtx", spliced_sparse)
    scipy.io.mmwrite(out_dir / f"unspliced_seed{seed}.mtx", unspliced_sparse)
    data.write_h5ad(out_dir / f"adata_seed{seed}.h5ad")
