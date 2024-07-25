import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import VELOVI
import datetime

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
savedata_folder = data_folder+"v2_pancreasINC/velovi/wopreprocess/"

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

########################################
########### split1
data_version = "split1"
print_message_with_time("#################### "+data_version+": Read data ")
adata = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_'+data_version+'_allgenes.h5ad') 

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("#################### "+data_version+": Preprocess data ")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### "+data_version+": Train model ")
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()

# save vae
print_message_with_time("#################### "+data_version+": Save vae ")
vae.save(savedata_folder+'vae_panINC_velovi_'+data_version+'_wopreprocess_v2.pt',overwrite=True)

# save original counts
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
adata.layers['spliced_original'] = S_mat[:,positions]
adata.layers['unspliced_original'] = U_mat[:,positions]

# write data
print_message_with_time("#################### "+data_version+": Save adata (final version) ")
adata.write(filename=savedata_folder+"adata_panINC_velovi_"+data_version+"_wopreprocess_v2.h5ad")
print_message_with_time("#################### "+data_version+": All done for "+data_version)

########################################
########### split1
data_version = "split2"
print_message_with_time("#################### "+data_version+": Read data ")
adata = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_'+data_version+'_allgenes.h5ad') 

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("#################### "+data_version+": Preprocess data ")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### "+data_version+": Train model ")
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()
# save vae
print_message_with_time("#################### "+data_version+": Save vae ")
vae.save(savedata_folder+'vae_panINC_velovi_'+data_version+'_wopreprocess_v2.pt',overwrite=True)
# save original counts
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
adata.layers['spliced_original'] = S_mat[:,positions]
adata.layers['unspliced_original'] = U_mat[:,positions]
# write data
print_message_with_time("#################### "+data_version+": Save adata (final version) ")
adata.write(filename=savedata_folder+"adata_panINC_velovi_"+data_version+"_wopreprocess_v2.h5ad")
print_message_with_time("#################### "+data_version+": All done for "+data_version)

########################################
########### total
data_version = "total"
print_message_with_time("#################### "+data_version+": Read data ")
adata = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_'+data_version+'_allgenes.h5ad') 

gene_names = adata.var.index.copy()
S_mat = adata.layers['spliced'].copy()
U_mat = adata.layers['unspliced'].copy()

print_message_with_time("#################### "+data_version+": Preprocess data ")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) # 30 in tutorial
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# train and apply model
print_message_with_time("#################### "+data_version+": Train model ")
VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()
# save vae
print_message_with_time("#################### "+data_version+": Save vae ")
vae.save(savedata_folder+'vae_panINC_velovi_'+data_version+'_wopreprocess_v2.pt',overwrite=True)
# save original counts
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}
positions = [positions_dict[gene] for gene in adata.var.index]
adata.layers['spliced_original'] = S_mat[:,positions]
adata.layers['unspliced_original'] = U_mat[:,positions]
# write data
print_message_with_time("#################### "+data_version+": Save adata (final version) ")
adata.write(filename=savedata_folder+"adata_panINC_velovi_"+data_version+"_wopreprocess_v2.h5ad")
print_message_with_time("#################### "+data_version+": All done for "+data_version)

##########################################################
##########################################################
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_velovi import *
method = 'velovi'
dataset_short = 'panINC'
dataset_long = 'pancreasINC'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/wopreprocess/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_v2.h5ad")
vae_split1 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_panINC_velovi_split1_wopreprocess_v2.pt', split1)

split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_v2.h5ad")
vae_split2 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_panINC_velovi_split2_wopreprocess_v2.pt', split2)

total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_v2.h5ad")
vae_total = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/wopreprocess/vae_panINC_velovi_total_wopreprocess_v2.pt', total)

raw = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
cell_index = np.array(np.where(raw.obs['clusters']!='Pre-endocrine')[0])
def create_adata_INC(S,U,adata_old):
    adata_new = ad.AnnData(X=S.astype(np.float32))
    adata_new.layers["spliced"] = S
    adata_new.layers["unspliced"] = U
    adata_new.uns = {}
    clusters_colors = dict(zip(adata_old.obs['clusters'].cat.categories,adata_old.uns['clusters_colors']))
    del clusters_colors['Pre-endocrine']
    adata_new.uns['clusters_colors'] = np.array(list(clusters_colors.values())).flatten().astype(object)
    adata_new.obs = adata_old.obs[adata_old.obs['clusters']!='Pre-endocrine']
    adata_new.obsm['X_pca'] = adata_old.obsm['X_pca'][cell_index,]
    adata_new.obsm['X_umap'] = adata_old.obsm['X_umap'][cell_index,]
    return adata_new

S_raw = raw.layers['spliced'][cell_index,:]
U_raw = raw.layers['unspliced'][cell_index,:]
raw = create_adata_INC(S=S_raw,U=U_raw,adata_old=raw)
total.uns['clusters_colors'] = split1.uns['clusters_colors'].copy()

## add velovi outputs to adata
add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)
## compute umap
compute_umap_pan(split1)
compute_umap_pan(split2)
compute_umap_pan(total)

split1.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split1_wopreprocess_outputAdded_v2.h5ad")
split2.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_split2_wopreprocess_outputAdded_v2.h5ad")
total.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/wopreprocess/adata_"+dataset_short+"_"+method+"_total_wopreprocess_outputAdded_v2.h5ad")
