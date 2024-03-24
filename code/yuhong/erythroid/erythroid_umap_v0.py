import scvelo as scv
import scanpy as sc
import bbknn

## total (the same for different seeds, so just use one of them)
### seed317
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317.h5ad')
scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_total, batch_key='sequencing.batch')
adata_total.X = adata_total.X.toarray()
bbknn.ridge_regression(adata_total, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_total)
bbknn.bbknn(adata_total, batch_key='sequencing.batch')
print("Batch correction done for total counts!")
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.velocity(adata_total)
scv.tl.velocity_graph(adata_total)
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_total_seed317_batchcorrected.png")

## split1
### seed317
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split1)
scv.pp.log1p(adata_split1)
scv.pp.moments(adata_split1, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_split1, batch_key='sequencing.batch')
adata_split1.X = adata_split1.X.toarray()
bbknn.ridge_regression(adata_split1, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split1)
bbknn.bbknn(adata_split1, batch_key='sequencing.batch')
print("Batch correction done for seed317 split1 counts!")
sc.pp.neighbors(adata_split1, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1)
scv.tl.velocity(adata_split1)
scv.tl.velocity_graph(adata_split1)
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed317_batchcorrected.png")
### seed320
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split1)
scv.pp.log1p(adata_split1)
scv.pp.moments(adata_split1, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_split1, batch_key='sequencing.batch')
adata_split1.X = adata_split1.X.toarray()
bbknn.ridge_regression(adata_split1, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split1)
bbknn.bbknn(adata_split1, batch_key='sequencing.batch')
print("Batch correction done for seed320 split1 counts!")
sc.pp.neighbors(adata_split1, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1)
scv.tl.velocity(adata_split1)
scv.tl.velocity_graph(adata_split1)
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed320_batchcorrected.png")


## split2
### seed317
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split2)
scv.pp.log1p(adata_split2)
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_split2, batch_key='sequencing.batch')
adata_split2.X = adata_split2.X.toarray()
bbknn.ridge_regression(adata_split2, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split2)
bbknn.bbknn(adata_split2, batch_key='sequencing.batch')
print("Batch correction done for seed317 split2 counts!")
sc.pp.neighbors(adata_split2, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2)
scv.tl.velocity(adata_split2)
scv.tl.velocity_graph(adata_split2)
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed317_batchcorrected.png")
### seed320
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split2)
scv.pp.log1p(adata_split2)
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)
### batch correction
bbknn.bbknn(adata_split2, batch_key='sequencing.batch')
adata_split2.X = adata_split2.X.toarray()
bbknn.ridge_regression(adata_split2, batch_key='sample', confounder_key='celltype')
sc.tl.pca(adata_split2)
bbknn.bbknn(adata_split2, batch_key='sequencing.batch')
print("Batch correction done for seed320 split2 counts!")
sc.pp.neighbors(adata_split2, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2)
scv.tl.velocity(adata_split2)
scv.tl.velocity_graph(adata_split2)
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed320_batchcorrected.png")

