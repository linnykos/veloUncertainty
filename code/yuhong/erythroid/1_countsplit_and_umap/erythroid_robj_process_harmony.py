import scvelo as scv
import scanpy as sc
import scanpy.external as sce

## total counts
print("read in total count data")
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317.h5ad')
adata_total.obs['sequencing.batch'] = adata_total.obs['sequencing.batch'].astype('category')
scv.pp.normalize_per_cell(adata_total)
scv.pp.log1p(adata_total)
scv.pp.moments(adata_total, n_pcs=30, n_neighbors=30)
### batch correction
sc.tl.pca(adata_total)
adata_total.obs['sequencing.batch'] = adata_total.obs['sequencing.batch'].astype('category')
sce.pp.harmony_integrate(adata_total, 'sequencing.batch', max_iter_harmony=20)
print("Batch correction done for total counts!")
sc.pp.neighbors(adata_total, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_total)
scv.tl.recover_dynamics(adata_total)
scv.tl.velocity(adata_total, mode="dynamical")
scv.tl.velocity_graph(adata_total)
### write h5ad object
adata_total.__dict__['_raw'].__dict__['_var'] = adata_total.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_total.raw)
#del(adata_total.var['_index']) 
adata_total.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317_harmony.h5ad")

## seed317 - split1
print("read in seed317 split1 count data")
adata_split1_317 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split1_317)
scv.pp.log1p(adata_split1_317)
scv.pp.moments(adata_split1_317, n_pcs=30, n_neighbors=30)
### batch correction
sc.tl.pca(adata_split1_317)
adata_split1_317.obs['sequencing.batch'] = adata_split1_317.obs['sequencing.batch'].astype('category')
sce.pp.harmony_integrate(adata_split1_317, 'sequencing.batch')
print("Batch correction done for seed317 split1 counts!")
sc.pp.neighbors(adata_split1_317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_317)
scv.tl.recover_dynamics(adata_split1_317)
scv.tl.velocity(adata_split1_317, mode="dynamical")
scv.tl.velocity_graph(adata_split1_317)
### write h5ad object
adata_split1_317.__dict__['_raw'].__dict__['_var'] = adata_split1_317.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split1_317.raw)
adata_split1_317.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317_harmony.h5ad")

## seed320 - split1
print("read in seed320 split1 count data")
adata_split1_320 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split1_320)
scv.pp.log1p(adata_split1_320)
scv.pp.moments(adata_split1_320, n_pcs=30, n_neighbors=30)
### batch correction
sc.tl.pca(adata_split1_320)
adata_split1_320.obs['sequencing.batch'] = adata_split1_320.obs['sequencing.batch'].astype('category')
sce.pp.harmony_integrate(adata_split1_320, 'sequencing.batch')
print("Batch correction done for seed320 split1 counts!")
sc.pp.neighbors(adata_split1_320, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split1_320)
scv.tl.recover_dynamics(adata_split1_320)
scv.tl.velocity(adata_split1_320, mode="dynamical")
scv.tl.velocity_graph(adata_split1_320)
### write h5ad object
adata_split1_320.__dict__['_raw'].__dict__['_var'] = adata_split1_320.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split1_320.raw)
adata_split1_320.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320_harmony.h5ad")

## seed317 - split2
print("read in seed317 split2 count data")
adata_split2_317 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317.h5ad")
scv.pp.normalize_per_cell(adata_split2_317)
scv.pp.log1p(adata_split2_317)
scv.pp.moments(adata_split2_317, n_pcs=30, n_neighbors=30)
### batch correction
sc.tl.pca(adata_split2_317)
adata_split2_317.obs['sequencing.batch'] = adata_split2_317.obs['sequencing.batch'].astype('category')
sce.pp.harmony_integrate(adata_split2_317, 'sequencing.batch')
print("Batch correction done for seed317 split2 counts!")
sc.pp.neighbors(adata_split2_317, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_317)
scv.tl.recover_dynamics(adata_split2_317)
scv.tl.velocity(adata_split2_317, mode="dynamical")
scv.tl.velocity_graph(adata_split2_317)
### write h5ad object
adata_split2_317.__dict__['_raw'].__dict__['_var'] = adata_split2_317.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split2_317.raw)
adata_split2_317.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317_harmony.h5ad")

## seed320 - split2
print("read in seed320 split2 count data")
adata_split2_320 = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320.h5ad")
scv.pp.normalize_per_cell(adata_split2_320)
scv.pp.log1p(adata_split2_320)
scv.pp.moments(adata_split2_320, n_pcs=30, n_neighbors=30)
### batch correction
sc.tl.pca(adata_split2_320)
adata_split2_320.obs['sequencing.batch'] = adata_split2_320.obs['sequencing.batch'].astype('category')
sce.pp.harmony_integrate(adata_split2_320, 'sequencing.batch')
print("Batch correction done for seed320 split2 counts!")
sc.pp.neighbors(adata_split2_320, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_split2_320)
scv.tl.recover_dynamics(adata_split2_320)
scv.tl.velocity(adata_split2_320, mode="dynamical")
scv.tl.velocity_graph(adata_split2_320)
### write h5ad object
adata_split2_320.__dict__['_raw'].__dict__['_var'] = adata_split2_320.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split2_320.raw)
adata_split2_320.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320_harmony.h5ad")






