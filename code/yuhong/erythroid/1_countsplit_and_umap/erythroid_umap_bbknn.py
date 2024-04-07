import scvelo as scv
import scanpy as sc

## total (the same for different seeds, so just use one of them)
## batch corrected
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317_bbknn.h5ad')
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_total_bbknn_celltype.png")
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_total_bbknn_seqbat.png")
print("**************** total counts done ****************")

## split1
### seed317
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317_bbknn.h5ad")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split1_seed317_bbknn_celltype.png")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split1_seed317_bbknn_seqbat.png")
print("**************** seed317 split1 done ****************")

### seed320
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320_bbknn.h5ad")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split1_seed320_bbknn_celltype.png")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split1_seed320_bbknn_seqbat.png")
print("**************** seed320 split1 done ****************")


## split2
### seed317
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317_bbknn.h5ad")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split2_seed317_bbknn_celltype.png")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split2_seed317_bbknn_seqbat.png")
print("**************** seed317 split2 done ****************")

### seed320
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320_bbknn.h5ad")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split2_seed320_bbknn_celltype.png")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/bbknn/scvelo_erythroid_split2_seed320_bbknn_seqbat.png")
print("**************** seed320 split2 done ****************")

