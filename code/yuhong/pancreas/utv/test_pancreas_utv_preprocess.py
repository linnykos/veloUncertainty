import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os

### did not run pp.neighbors
# the below script uses the environment: "utvClone"

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = False
velo_config.IROOT = None
velo_config.FIT_OPTION = '2'
velo_config.ASSIGN_POS_U = True

os.environ["TF_USE_LEGACY_KERAS"]="1"

label='clusters'
adata = scv.datasets.pancreas()

adata_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad"
adata_res = utv.run_model(adata_path, label, config_file=velo_config)
sc.tl.umap(adata_res)
adata.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_utv/pancreas_utv_preprocess.h5ad")
scv.pl.velocity_embedding_stream(adata_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/test_unitvelo_preprocess.png")

