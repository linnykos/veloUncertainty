# https://www.sciencedirect.com/science/article/pii/S2589004222018260

target_genes = ["1110034G24Rik", "Abcg2", "Arl4a", "Cdc42ep3", "Cpeb4", "Gda", "Ghitm", "Hace1", "Hpf1", "Ifit1bl1", "Kel", "Lrrc1", "Mllt3", "Ncoa7", "Pkhd1l1", "Pla2g16", "Rhag", "Scai", "Sec61g", "Smc2os", "Snca", "Trim10", "Tspan8", "Zfand6", "Zfp760"]
np.intersect1d(total.var.index, target_genes)

sc.tl.rank_genes_groups(total, 'celltype', method='t-test')


adata = total[:,target_genes].copy()

scv.pp.filter_and_normalize(adata, min_shared_counts=20)
sc.tl.pca(adata, svd_solver="arpack")


scv.pp.moments(adata, n_pcs=5, n_neighbors=30)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=5)
sc.tl.umap(adata)
scv.tl.recover_dynamics(adata,n_jobs=1)
scv.tl.velocity(adata, mode="dynamical", min_r2=0, min_likelihood=0) # https://scvelo.readthedocs.io/en/stable/scvelo.tl.velocity.html
scv.tl.velocity_graph(adata,n_jobs=1)