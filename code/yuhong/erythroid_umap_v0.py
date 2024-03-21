import scvelo as scv
import scanpy as sc

## total
### seed317
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317.h5ad')
scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_total, svd_solver="arpack")
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.velocity(adata_total)
scv.tl.velocity_graph(adata_total)
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo_erythroid_total_seed317.png")
### seed320
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed320.h5ad')
scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_total, svd_solver="arpack")
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.velocity(adata_total)
scv.tl.velocity_graph(adata_total)
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo_erythroid_total_seed320.png")


## split1
### seed317
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split1)
scv.pp.log1p(adata_split1)
scv.pp.moments(adata_split1, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1, svd_solver="arpack")
sc.pp.neighbors(adata_split1, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1)
scv.tl.velocity(adata_split1)
scv.tl.velocity_graph(adata_split1)
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo_erythroid_split1_seed317.png")
### seed320
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split1)
scv.pp.log1p(adata_split1)
scv.pp.moments(adata_split1, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1, svd_solver="arpack")
sc.pp.neighbors(adata_split1, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1)
scv.tl.velocity(adata_split1)
scv.tl.velocity_graph(adata_split1)
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo_erythroid_split1_seed320.png")


## split2
### seed317
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split2)
scv.pp.log1p(adata_split2)
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2, svd_solver="arpack")
sc.pp.neighbors(adata_split2, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2)
scv.tl.velocity(adata_split2)
scv.tl.velocity_graph(adata_split2)
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo_erythroid_split2_seed317.png")
### seed320
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split2)
scv.pp.log1p(adata_split2)
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2, svd_solver="arpack")
sc.pp.neighbors(adata_split2, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2)
scv.tl.velocity(adata_split2)
scv.tl.velocity_graph(adata_split2)
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/scvelo_erythroid_split2_seed320.png")

