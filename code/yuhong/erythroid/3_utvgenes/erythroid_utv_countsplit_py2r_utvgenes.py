import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
from scipy.sparse import csr_matrix

### did not run pp.neighbors
# the below script uses the environment: "utvClone"

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1
#velo_config.ASSIGN_POS_U = True

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

spliced = adata.layers['spliced'].copy() # shape=(9815, 53801)
unspliced = adata.layers['unspliced'].copy()
gene_names = adata.var['Accession'].copy()

adata.X = csr_matrix(adata.X)
positions_dict = {gene: pos for pos, gene in enumerate(gene_names.index)}

positions = [positions_dict[gene] for gene in adata.var['Accession'].index]

spliced_subset = spliced[:,positions]
unspliced_subset = unspliced[:,positions]
adata.layers['spliced_original'] = spliced_subset
adata.layers['unspliced_original'] = unspliced_subset

path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup6_eryutv_utvgenes/tmp_ery_utvgenes.h5ad"
adata.write_h5ad(filename=path)
adata_res = utv.run_model(adata, label, config_file=velo_config)
adata_res.write_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utv_preprocess.h5ad')

scv.pl.velocity_embedding_stream(adata_res,basis="umap",color="celltype",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/unitvelo_preprocess.png")
scv.pl.velocity_embedding_stream(adata_res,color=adata.uns['label'],dpi=100,
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/unitvelo_preprocess_tut.png")

