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
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_erythroid/utv/" 

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

###########
print_message_with_time("#################### Read data ")
adata_split1 = sc.read_h5ad(data_folder+'v2_erythroid/seed317_split1_allgenes.h5ad')

gene_names = adata_split1.var.index.copy()
S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()

### fit model
print_message_with_time("#################### Fit model")
adata_split1.write_h5ad(filename="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup16_v2_utv/tmp_v2_ery_split1.h5ad")
adata_split1 = utv.run_model("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup16_v2_utv/tmp_v2_ery_split1.h5ad", 'celltype', config_file=velo_config)

print_message_with_time("#################### Write data (intermediate version)")
adata_split1.write_h5ad(data_folder+'v2_erythroid/utv/backup/adata_ery_utv_split1.h5ad')
### 
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions_split1 = [positions_dict[gene] for gene in adata_split1.var.index]

adata_split1.layers['spliced_original'] = S_mat_split1[:,positions_split1] # looks like i did not do this actually
adata_split1.layers['unspliced_original'] = U_mat_split1[:,positions_split1]

print_message_with_time("#################### Write data (final version)")
adata_split1.write_h5ad(data_folder+'v2_erythroid/utv/adata_ery_utv_split1_v2.h5ad')

print_message_with_time("#################### All done for split1")



