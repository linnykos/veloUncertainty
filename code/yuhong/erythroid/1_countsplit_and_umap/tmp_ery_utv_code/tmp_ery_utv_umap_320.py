import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os

# the below script uses the environment: "utvClone"

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1

os.environ["TF_USE_LEGACY_KERAS"]="1"
# https://github.com/tensorflow/tensorflow/issues/62337
# To solve the following error: 
## ImportError: `keras.optimizers.legacy` is not supported in Keras 3. 
## When using `tf.keras`, to continue using a `tf.keras.optimizers.legacy` optimizer, 
## you can install the `tf_keras` package (Keras 2) and set the environment variable `TF_USE_LEGACY_KERAS=True` 
## to configure TensorFlow to use `tf_keras` when accessing `tf.keras`.
# Tried this but did not work:
## velo_config.TF_USE_LEGACY_KERAS=True

label='celltype'
adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad")
#adata = scv.datasets.gastrulation_erythroid()
scv.pp.filter_genes(adata, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) 
### Extracted 2000 highly variable genes.
scv.pp.log1p(adata)

true_indices = adata.var['highly_variable'][adata.var['highly_variable'] == True].index.tolist()
velo_config.VGENES = true_indices


split1_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed320_split1_seurat.h5ad')
split1_seed320.var['highly_variable'] = adata.var['highly_variable'].copy()
split1_seed320.__dict__['_raw'].__dict__['_var'] = split1_seed320.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split1_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup3_eryutv/out/utv_erythroid_split1_seurat_seed320.h5ad"
split1_seed320.write_h5ad(filename=split1_path)
split1_seed320_res = utv.run_model(split1_path, label, config_file=velo_config)
split1_seed320_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed320_split1.h5ad')
sc.tl.umap(split1_seed320_res)
scv.pl.velocity_embedding_stream(split1_seed320_res,basis="umap",color=label,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_uncorrected_seed320_split1.png")
split1_seed320_res.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(split1_seed320_res,basis="umap",color=split1_seed320_res.uns['label'],
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_preumap_seed320_split1.png")
print("************* write split1_seed320_res to '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed320_split1.h5ad' *******************")


split2_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/erythroid_seed320_split2_seurat.h5ad')
split2_seed320.var['highly_variable'] = adata.var['highly_variable'].copy()
split2_seed320.__dict__['_raw'].__dict__['_var'] = split2_seed320.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split2_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup3_eryutv/out/utv_erythroid_split2_seurat_seed320.h5ad"
split2_seed320.write_h5ad(filename=split2_path)
split2_seed320_res = utv.run_model(split2_path, label, config_file=velo_config)
split2_seed320_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed320_split2.h5ad')
sc.tl.umap(split2_seed320_res)
scv.pl.velocity_embedding_stream(split2_seed320_res,basis="umap",color=label,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_uncorrected_seed320_split2.png")
split2_seed320_res.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(split2_seed320_res,basis="umap",color=split2_seed320_res.uns['label'],
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo/unitvelo_preumap_seed320_split2.png")
print("************* write split2_seed320_res to '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_utv/erythroid_utv_seed320_split2.h5ad' *******************")
