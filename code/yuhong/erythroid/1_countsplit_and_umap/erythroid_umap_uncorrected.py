import scvelo as scv
import scanpy as sc

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

adata = scv.datasets.gastrulation_erythroid()

scv.pp.filter_genes(adata, min_shared_counts=20)
### Filtered out 47456 genes that are detected 20 counts (shared).
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
### Extracted 2000 highly variable genes.
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_preprocess_uncorrected_celltype.png")
scv.pl.velocity_embedding_stream(adata, basis='umap', color='sequencing.batch',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_preprocess_uncorrected_seqbatch.png")
print("####### Uncorrected UMAP for raw data (2000 genes filtered) plotted #######")

## total (the same for different seeds, so just use one of them)
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
scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_preprocess_uncorrected_celltype.png")
scv.pl.velocity_embedding_stream(adata, basis='umap', color='sequencing.batch',save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_preprocess_uncorrected_seqbatch.png")
print("####### Uncorrected UMAP for total count data (converted through R) plotted #######")

## split1
### seed317
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split1)
scv.pp.log1p(adata_split1)
scv.pp.moments(adata_split1, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1)
sc.pp.neighbors(adata_split1, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1)
scv.tl.recover_dynamics(adata_split1)
scv.tl.velocity(adata_split1, mode="dynamical")
scv.tl.velocity_graph(adata_split1)
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed317_uncorrected_celltype.png")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed317_uncorrected_seqbat.png")

### seed320
adata_split1 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split1)
scv.pp.log1p(adata_split1)
scv.pp.moments(adata_split1, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split1)
sc.pp.neighbors(adata_split1, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1)
scv.tl.recover_dynamics(adata_split1)
scv.tl.velocity(adata_split1, mode="dynamical")
scv.tl.velocity_graph(adata_split1)
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed320_uncorrected_celltype.png")
scv.pl.velocity_embedding_stream(adata_split1, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split1_seed320_uncorrected_seqbat.png")


## split2
### seed317
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split2)
scv.pp.log1p(adata_split2)
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2)
sc.pp.neighbors(adata_split2, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2)
scv.tl.recover_dynamics(adata_split2)
scv.tl.velocity(adata_split2, mode="dynamical")
scv.tl.velocity_graph(adata_split2)
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed317_uncorrected_celltype.png")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed317_uncorrected_seqbat.png")

### seed320
adata_split2 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split2)
scv.pp.log1p(adata_split2)
scv.pp.moments(adata_split2, n_pcs=30, n_neighbors=30)
sc.tl.pca(adata_split2)
sc.pp.neighbors(adata_split2, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2)
scv.tl.recover_dynamics(adata_split2)
scv.tl.velocity(adata_split2, mode="dynamical")
scv.tl.velocity_graph(adata_split2)
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed320_uncorrected_celltype.png")
scv.pl.velocity_embedding_stream(adata_split2, basis='umap',color="sequencing.batch",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo_erythroid_split2_seed320_uncorrected_seqbat.png")





