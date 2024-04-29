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
adata = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utv_preprocess.h5ad')

split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/erythroid_seed317_split1_seurat.h5ad')
split1_seed317.var['highly_variable'] = adata.var['highly_variable'].copy()
split1_seed317.__dict__['_raw'].__dict__['_var'] = split1_seed317.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split1_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup6_eryutv_utvgenes/utv_erythroid_split1_seurat_seed317.h5ad"
split1_seed317.write_h5ad(filename=split1_path)
split1_seed317_res = utv.run_model(split1_path, label, config_file=velo_config)
split1_seed317_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utvgenes_seed317_split1.h5ad')
# using original adata umap
split1_seed317_res2 = split1_seed317_res.copy()
split1_seed317_res2.obsm['X_umap'] = adata.obsm['X_umap'].copy()
scv.pl.velocity_embedding_stream(split1_seed317_res2,basis="umap",color=split1_seed317_res2.uns['label'],
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_preumap_317s1.png")
# uncorrected recovered umap
sc.tl.umap(split1_seed317_res)
scv.pl.velocity_embedding_stream(split1_seed317_res,basis="umap",color=label,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_uncorrected_317s1.png")
### batch correction
import bbknn
bbknn.bbknn(split1_seed317_res, batch_key='sequencing.batch')
split1_seed317_res.X = split1_seed317_res.X.toarray()
bbknn.ridge_regression(split1_seed317_res, batch_key='sample', confounder_key='celltype')
sc.tl.pca(split1_seed317_res)
bbknn.bbknn(split1_seed317_res, batch_key='sequencing.batch')
sc.pp.neighbors(split1_seed317_res, n_neighbors=10, n_pcs=40)
sc.tl.umap(split1_seed317_res)
scv.pl.velocity_embedding_stream(split1_seed317_res,basis="umap",color=label,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/ery_utvgenes_bbknn_317s1.png")
print("************* write split1_seed317_res to '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/' *******************")
