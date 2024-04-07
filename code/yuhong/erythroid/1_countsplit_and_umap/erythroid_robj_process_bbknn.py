import scvelo as scv
import scanpy as sc
import bbknn

### total counts
print("********************* Read in total count data **********************")
adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317_Rbbknn.h5ad')
adata_total.obs['sequencing.batch'] = adata_total.obs['sequencing.batch'].astype('category')
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
scv.tl.recover_dynamics(adata_total)
scv.tl.velocity(adata_total, mode="dynamical")
scv.tl.velocity_graph(adata_total)
### write h5ad object
adata_total.__dict__['_raw'].__dict__['_var'] = adata_total.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_total.raw)
#del(adata_total.var['_index']) 
adata_total.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_total_seurat_seed317_bbknn.h5ad")
adata_total = ''
print("********************* Wrote total count data **********************")


## seed317, split1
print("********************* Read in seed317 split1 data **********************")
adata_split1 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317_Rbbknn.h5ad')
adata_split1.obs['sequencing.batch'] = adata_split1.obs['sequencing.batch'].astype('category')
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
scv.tl.recover_dynamics(adata_split1)
scv.tl.velocity(adata_split1, mode="dynamical")
scv.tl.velocity_graph(adata_split1)
### write h5ad object
adata_split1.__dict__['_raw'].__dict__['_var'] = adata_split1.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split1.raw)
#del(adata_split1.var['_index']) 
adata_split1.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed317_bbknn.h5ad")
adata_split1=''
print("********************* Wrote seed317 split1 data **********************")



## seed317, split2
print("********************* Read in seed317 split2 data **********************")
adata_split2 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317_Rbbknn.h5ad')
adata_split2.obs['sequencing.batch'] = adata_split2.obs['sequencing.batch'].astype('category')
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
scv.tl.recover_dynamics(adata_split2)
scv.tl.velocity(adata_split2, mode="dynamical")
scv.tl.velocity_graph(adata_split2)
### write h5ad object
adata_split2.__dict__['_raw'].__dict__['_var'] = adata_split2.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split2.raw)
#del(adata_total.var['_index']) 
adata_split2.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed317_bbknn.h5ad")
adata_split2=''
print("********************* Wrote seed317 split2 data **********************")


## seed320, split1
print("********************* Read in seed317 split1 data **********************")
adata_split1 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320_Rbbknn.h5ad')
adata_split1.obs['sequencing.batch'] = adata_split1.obs['sequencing.batch'].astype('category')
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
scv.tl.recover_dynamics(adata_split1)
scv.tl.velocity(adata_split1, mode="dynamical")
scv.tl.velocity_graph(adata_split1)
### write h5ad object
adata_split1.__dict__['_raw'].__dict__['_var'] = adata_split1.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split1.raw)
#del(adata_split1.var['_index']) 
adata_split1.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split1_seurat_seed320_bbknn.h5ad")
adata_split1=''
print("********************* Wrote seed320 split1 data **********************")


## seed317, split2
print("********************* Read in seed317 split2 data **********************")
adata_split2 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320_Rbbknn.h5ad')
adata_split2.obs['sequencing.batch'] = adata_split2.obs['sequencing.batch'].astype('category')
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
scv.tl.recover_dynamics(adata_split2)
scv.tl.velocity(adata_split2, mode="dynamical")
scv.tl.velocity_graph(adata_split2)
### write h5ad object
adata_split2.__dict__['_raw'].__dict__['_var'] = adata_split2.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata_split2.raw)
#del(adata_total.var['_index']) 
adata_split2.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/scvelo_erythroid_split2_seurat_seed320_bbknn.h5ad")
adata_split2=''
print("********************* Wrote seed320 split2 data **********************")

