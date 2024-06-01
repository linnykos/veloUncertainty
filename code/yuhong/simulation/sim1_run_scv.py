import scvelo as scv
import scanpy as sc

data1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim1.h5ad")
scv.pp.normalize_per_cell(data1,enforce=True)
scv.pp.log1p(data1)
scv.pp.moments(data1, n_pcs=30, n_neighbors=30)
sc.tl.pca(data1, svd_solver="arpack")
sc.pp.neighbors(data1, n_neighbors=10, n_pcs=40)
sc.tl.umap(data1)
scv.tl.recover_dynamics(data1)
scv.tl.velocity(data1, mode="dynamical")
scv.tl.velocity_graph(data1)
scv.pl.velocity_embedding_stream(data1, basis='pca', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1_pca_vcompute.png")
scv.pl.velocity_embedding_stream(data1, basis='umap', save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1_umap_vcompute.png")
scv.tl.velocity_confidence(data1)
scv.pl.scatter(data1, c='velocity_confidence', basis='pca',cmap='coolwarm', perc=[5, 95],
               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1_pca_vcompute_scvconf.png")
scv.pl.scatter(data1, c='velocity_confidence', basis='umap',cmap='coolwarm', perc=[5, 95],
               save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1_umap_vcompute_scvconf.png")


# alpha = 5, beta = .3, gamma = .5
data1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim1.h5ad")
scv.pp.normalize_per_cell(data1,enforce=True)
scv.pp.log1p(data1)
scv.pp.moments(data1, n_pcs=30, n_neighbors=30)
sc.tl.pca(data1, svd_solver="arpack")
sc.pp.neighbors(data1, n_neighbors=10, n_pcs=40)
sc.tl.umap(data1,alpha=5,gamma=.5)
scv.tl.recover_dynamics(data1,t_max=25)
scv.tl.velocity_graph(data1,vkey="true_velocity")
scv.pl.velocity_embedding_stream(data1, basis='pca',vkey="true_velocity", save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1_pca_vtrue.png")
scv.pl.velocity_embedding_stream(data1, basis='umap',vkey="true_velocity", save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim1_umap_vtrue.png")

