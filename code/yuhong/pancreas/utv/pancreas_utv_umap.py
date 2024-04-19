import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
from sklearn.metrics.pairwise import cosine_similarity

### did not run pp.neighbors
# the below script uses the environment: "utvClone"

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = False
velo_config.IROOT = None
velo_config.FIT_OPTION = '2'
velo_config.ASSIGN_POS_U = True

os.environ["TF_USE_LEGACY_KERAS"]="1"
# https://github.com/tensorflow/tensorflow/issues/62337
# To solve the following error: 
## ImportError: `keras.optimizers.legacy` is not supported in Keras 3. 
## When using `tf.keras`, to continue using a `tf.keras.optimizers.legacy` optimizer, 
## you can install the `tf_keras` package (Keras 2) and set the environment variable `TF_USE_LEGACY_KERAS=True` 
## to configure TensorFlow to use `tf_keras` when accessing `tf.keras`.
# Tried this but did not work:
## velo_config.TF_USE_LEGACY_KERAS=True

label='clusters'
adata = scv.datasets.pancreas()
scv.pp.filter_genes(adata, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) 

true_indices = adata.var['highly_variable_genes'][adata.var['highly_variable_genes'] == 'True'].index.tolist()
velo_config.VGENES = true_indices

adata_total = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed317_total_seurat.h5ad')
adata_total.var['highly_variable'] = adata.var['highly_variable'].copy()
adata_total.__dict__['_raw'].__dict__['_var'] = adata_total.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# To solve the error: ValueError: '_index' is a reserved name for dataframe columns.
# ref: https://github.com/theislab/scvelo/issues/255
total_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup2_utv/tmp/scvelo_pancreas_total_seurat_seed317.h5ad"
adata_total.write_h5ad(filename=total_path)
adata_total_res = utv.run_model(total_path, label, config_file=velo_config)
sc.tl.umap(adata_total_res)
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_total.png")
print("total counts done!")

split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed317_split1_seurat.h5ad')
split1_seed317.var['highly_variable'] = adata.var['highly_variable'].copy()
split1_seed317.__dict__['_raw'].__dict__['_var'] = split1_seed317.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split1_seed317_path = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup2_utv/tmp/scvelo_pancreas_split1_seurat_seed317.h5ad'
split1_seed317.write_h5ad(filename=split1_seed317_path)
split1_seed317_res = utv.run_model(split1_seed317_path, label, config_file=velo_config)
sc.tl.umap(split1_seed317_res)
scv.pl.velocity_embedding_stream(split1_seed317_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_split1.png")
print("seed317 split1 done!")

split2_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed317_split2_seurat.h5ad')
split2_seed317.var['highly_variable'] = adata.var['highly_variable'].copy()
split2_seed317.__dict__['_raw'].__dict__['_var'] = split2_seed317.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split2_seed317_path = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup2_utv/tmp/scvelo_pancreas_split2_seurat_seed317.h5ad'
split2_seed317.write_h5ad(filename=split2_seed317_path)
split2_seed317_res = utv.run_model(split2_seed317_path, label, config_file=velo_config)
sc.tl.umap(split2_seed317_res)
scv.pl.velocity_embedding_stream(split2_seed317_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_split2.png")
print("seed317 split2 done!")

split1_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed320_split1_seurat.h5ad')
split1_seed320.var['highly_variable'] = adata.var['highly_variable'].copy()
split1_seed320.__dict__['_raw'].__dict__['_var'] = split1_seed320.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split1_seed320_path = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup2_utv/tmp/scvelo_pancreas_split1_seurat_seed320.h5ad'
split1_seed320.write_h5ad(filename=split1_seed320_path)
split1_seed320_res = utv.run_model(split1_seed320_path, label, config_file=velo_config)
sc.tl.umap(split1_seed320_res)
scv.pl.velocity_embedding_stream(split1_seed320_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed320_split1.png")
print("seed320 split1 done!")

split2_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed320_split2_seurat.h5ad')
split2_seed320.var['highly_variable'] = adata.var['highly_variable'].copy()
split2_seed320.__dict__['_raw'].__dict__['_var'] = split2_seed320.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split2_seed320_path = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup2_utv/tmp/scvelo_pancreas_split2_seurat_seed320.h5ad'
split2_seed320.write_h5ad(filename=split2_seed320_path)
split2_seed320_res = utv.run_model(split2_seed320_path, label, config_file=velo_config)
sc.tl.umap(split2_seed320_res)
scv.pl.velocity_embedding_stream(split2_seed320_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed320_split2.png")
print("seed320 split2 done!")

adata_total_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_total.h5ad')
split1_seed317_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed317_split1.h5ad')
split2_seed317_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed317_split2.h5ad')
split1_seed320_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed320_split1.h5ad')
split2_seed320_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_seed320_split2.h5ad')

