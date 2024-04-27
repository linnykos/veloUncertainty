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
velo_config.FIT_OPTION = '1'
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
adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utv_preprocess.h5ad")

split1_seed320 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split1_seurat.h5ad')
split1_seed320.var['highly_variable'] = adata.var['highly_variable'].copy()
split1_seed320.__dict__['_raw'].__dict__['_var'] = split1_seed320.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split1_seed320_path = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup5_panutv_utvgenes/tmp_pan_utvgenes_split1_seurat_seed320.h5ad'
split1_seed320.write_h5ad(filename=split1_seed320_path)
split1_seed320_res = utv.run_model(split1_seed320_path, label, config_file=velo_config)
sc.tl.umap(split1_seed320_res)
scv.pl.velocity_embedding_stream(split1_seed320_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/utvgenes_seed320s1.png")
print("seed320 split1 done!")

split1_seed320_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split1.h5ad')

