import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
import datetime

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

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

###########
data_version = 'split1'
print_message_with_time("#################### Read data ")
adata = sc.read_h5ad(data_folder+'v2_larry/larry_'+data_version+'_allgenes.h5ad')

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

### fit model
print_message_with_time("#################### Fit model")
adata = utv.run_model(data_folder+'v2_larry/larry_'+data_version+'_allgenes.h5ad', 'state_info', config_file=velo_config)

### 
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]

adata.layers['spliced_original'] = S_mat[:,positions] # looks like i did not do this actually
adata.layers['unspliced_original'] = U_mat[:,positions]

print_message_with_time("#################### Write data (final version)")
adata.write_h5ad(data_folder+'v2_larry/utv/adata_larry_utv_'+data_version+'_v2.h5ad')

print_message_with_time("#################### All done for split1")




