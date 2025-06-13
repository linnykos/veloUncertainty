import numpy as np
import scipy.sparse as sp
import anndata as ad
from typing import Optional, Tuple
import cellrank as cr


# ---------------------------------------------------------------------
# 1. Chung directed Laplacian (unchanged, but with a short docstring)
# ---------------------------------------------------------------------
def chung_laplacian(
    A_csr: sp.csr_matrix,
    tol: float = 1e-12,
    max_iter: int = 10_000,
) -> Tuple[sp.csr_matrix, np.ndarray]:
    """
    Compute Chung’s directed Laplacian L and the stationary distribution π
    for a (possibly non-normalised) transition matrix A.

    Parameters
    ----------
    A_csr
        Square sparse matrix (n × n) in CSR format whose rows sum to ≥ 0.
    tol
        1-norm convergence tolerance for power iteration on π.
    max_iter
        Maximum number of power-iteration steps.

    Returns
    -------
    L
        Chung directed Laplacian (symmetric, PSD, sparse CSR).
    pi
        Length-n stationary distribution (π ≥ 0, πᵀ 1 = 1).
    """
    # row-normalise to a stochastic matrix P
    d_out = np.asarray(A_csr.sum(axis=1)).ravel()
    P = sp.diags(1.0 / np.maximum(d_out, 1e-32)) @ A_csr

    # power iteration for πᵀ P = πᵀ
    pi = np.ones(P.shape[0]) / P.shape[0]
    for _ in range(max_iter):
        new_pi = pi @ P
        if np.linalg.norm(new_pi - pi, 1) < tol:
            break
        pi = new_pi

    S = sp.diags(np.sqrt(pi)) @ P @ sp.diags(1.0 / np.sqrt(pi))
    H = 0.5 * (S + S.T)
    L = sp.eye(P.shape[0], format="csr") - H
    return L, pi


# ---------------------------------------------------------------------
# 2. Align / subset data and graph
# ---------------------------------------------------------------------
def prepare_adata_laplacian(
    adata_kernel: ad.AnnData,
    adata_full: ad.AnnData,
    subset_key: Optional[str] = None,
    subset_value: Optional[Union[str, list[str], set[str]]] = None,
) -> Tuple[ad.AnnData, sp.csr_matrix]:
    """
    Align cell order, optionally subset, and return the trimmed AnnData
    plus its corresponding velocity graph.

    Parameters
    ----------
    adata_kernel
        The AnnData used to build the CellRank/VelocityKernel.
    adata_full
        The AnnData holding all genes to be scored.
    subset_key, subset_value
        If provided, keep only cells where
        `adata_kernel.obs[subset_key] ∈ subset_value`.

    Returns
    -------
    adata_trim
        `adata_full` restricted to the chosen cells, with identical order
        as the returned graph.
    vel_trim
        Velocity graph restricted to (and re-ordered with) the same cells.
    """
    # 1) intersect cells and bring both AnnData objects into the same order
    shared = np.intersect1d(adata_kernel.obs_names, adata_full.obs_names)
    adata_kernel = adata_kernel[shared].copy()
    adata_full   = adata_full[shared].copy()
    
    # 2) Compute the velocity graph
    vk = cr.kernels.VelocityKernel(adata_kernel)
    vk.compute_transition_matrix()
    velocity_graph = vk.transition_matrix

    # 2) optional sub-selection
    if subset_key is not None and subset_value is not None:
        if not isinstance(subset_value, (list, set, tuple)):
            subset_value = [subset_value]
        mask = adata_kernel.obs[subset_key].isin(subset_value)
        adata_kernel = adata_kernel[mask].copy()
        adata_full = adata_full[mask].copy()
        row_idx = np.where(mask)[0]          # indices to keep
        velocity_graph = velocity_graph[row_idx, :][:, row_idx]

    return adata_full, velocity_graph

# ---------------------------------------------------------------------
# 3. Laplacian-score calculation for every gene
# ---------------------------------------------------------------------
def compute_laplacian(
    adata: ad.AnnData,
    velocity_graph: sp.csr_matrix,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Compute Chung-Laplacian scores for **every** gene in `adata`.

    Returns
    -------
    scores
        1-D array (len = adata.n_vars); NaN for genes with zero variance.
    """
    L, pi = chung_laplacian(velocity_graph)
    n_genes = adata.n_vars
    scores = np.empty(n_genes, dtype=float)

    for g in range(n_genes):
        # expression vector x (shape n_cells,)
        vec = adata[:, g].X
        x = vec.toarray().ravel() if sp.issparse(vec) else np.asarray(vec).ravel()

        # numerator  xᵀ L x
        numer = x.dot(L.dot(x))

        # variance under π
        nz = x != 0
        mean_pi = np.dot(pi[nz], x[nz])
        var_pi = (
            np.dot(pi[nz], (x[nz] - mean_pi) ** 2)
            + (1.0 - pi[nz].sum()) * mean_pi**2
        )

        scores[g] = np.nan if var_pi < eps else numer / var_pi

    return scores


# ---------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------
if __name__ == "__main__":
    import cellrank as cr
    import scanpy as sc

    # --- 1. load raw objects (paths omitted for brevity) -------------
    adata_kernel = ad.read_h5ad("path/to/adata_velo.h5ad")
    adata_full = ad.read_h5ad("path/to/adata_full.h5ad")

    # --- 2. build velocity graph ------------------------------------
    vk = cr.kernels.VelocityKernel(adata_kernel)
    vk.compute_transition_matrix()
    v_graph = vk.transition_matrix

    # --- 3. align / subset ------------------------------------------
    adata_sub, v_graph_sub = prepare_adata_laplacian(
        adata_kernel,
        adata_full,
        v_graph,
        subset_key="cluster_name",
        subset_value="excitatory neuron 4",
    )

    # --- 4. Laplacian scores ----------------------------------------
    lap_scores = compute_laplacian(adata_sub, v_graph_sub)
    print("Done – scores shape:", lap_scores.shape)
