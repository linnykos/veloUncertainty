import scvelo as scv
import scanpy as sc
import scanpy.external as sce

## total (the same for different seeds, so just use one of them)
### uncorrected
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317.h5ad')
scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_total)
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.recover_dynamics(adata_total)
scv.tl.velocity(adata_total, mode="dynamical")
scv.tl.velocity_graph(adata_total)
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_total_uncorrected_harmony_celltype.png")
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_total_uncorrected_harmony_seqbat.png")
print("################### Uncorrected UMAP for total counts plotted! ###################")

### batch corrected
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317_harmony.h5ad')
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_total_harmony_celltype.png")
scv.pl.velocity_embedding_stream(adata_total, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_total_harmony_seqbat.png")
print("################### Corrected UMAP for total counts plotted! ###################")


## split1
### seed317
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317_harmony.h5ad")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed317_harmony_celltype.png")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed317_harmony_seqbat.png")
print("################### Corrected UMAP for seed317 split1 counts plotted! ###################")


### seed320
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320_harmony.h5ad")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed320_harmony_celltype.png")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed320_harmony_seqbat.png")
print("################### Corrected UMAP for seed320 split1 counts plotted! ###################")


## split2
### seed317
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317_harmony.h5ad")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed317_harmony_celltype.png")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed317_harmony_seqbat.png")
print("################### Corrected UMAP for seed317 split2 counts plotted! ###################")


### seed320
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320_harmony.h5ad")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed320_harmony_celltype.png")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed320_harmony_seqbat.png")
print("################### Corrected UMAP for seed320 split2 counts plotted! ###################")

