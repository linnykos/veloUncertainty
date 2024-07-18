import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
from scipy.sparse import csr_matrix
import datetime

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

dataset_long = "pancreas"
dataset_short = "pan"
method = "utv"
celltype_label = "clusters"
data_version = "total"

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
#fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/" 

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

print_message_with_time("#################### Read data ")
adata = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad") # 3696 × 27998

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

### fit model
print_message_with_time("#################### Fit model")
adata.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup16_v2_utv/tmp_v2_"+dataset_short+"_total.h5ad")
adata = utv.run_model("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup16_v2_utv/tmp_v2_"+dataset_short+"_total.h5ad", celltype_label, config_file=velo_config)

### 
print_message_with_time("#################### Write original counts")
positions = [positions_dict[gene] for gene in adata.var.index]
adata.layers['spliced_original'] = S_mat[:,positions] 
adata.layers['unspliced_original'] = U_mat[:,positions]

print_message_with_time("#################### Write data")
adata.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+data_version+'_v2.h5ad')

print_message_with_time("####################All done for "+data_version)



