import scvelo as scv
import scanpy as sc

adata_pan_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/scvelo_pancreas_total_seurat.h5ad')

scv.pp.normalize_per_cell(adata_pan_total)
scv.pp.log1p(adata_pan_total)
scv.pp.moments(adata_pan_total, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_pan_total, svd_solver="arpack")
sc.pp.neighbors(adata_pan_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_pan_total)
scv.tl.velocity(adata_pan_total)
scv.tl.velocity_graph(adata_pan_total)
scv.pl.velocity_embedding_stream(adata_pan_total, basis='umap',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/pancreas_total_1.png")

adata_pan_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/scvelo_pancreas_split1_seurat.h5ad")
scv.pp.normalize_per_cell(adata_pan_split1)
scv.pp.log1p(adata_pan_split1)
scv.pp.moments(adata_pan_split1, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_pan_split1, svd_solver="arpack")
sc.pp.neighbors(adata_pan_split1, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_pan_split1)
scv.tl.velocity(adata_pan_split1)
scv.tl.velocity_graph(adata_pan_split1)
scv.pl.velocity_embedding_stream(adata_pan_split1, basis='umap',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/adata_pan_split1_1.png")


adata_pan_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/scvelo_pancreas_split2_seurat.h5ad")
scv.pp.normalize_per_cell(adata_pan_split2)
scv.pp.log1p(adata_pan_split2)
scv.pp.moments(adata_pan_split2, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_pan_split2, svd_solver="arpack")
sc.pp.neighbors(adata_pan_split2, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_pan_split2)
scv.tl.velocity(adata_pan_split2)
scv.tl.velocity_graph(adata_pan_split2)
scv.pl.velocity_embedding_stream(adata_pan_split2, basis='umap',
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/adata_pan_split2_1.png")

