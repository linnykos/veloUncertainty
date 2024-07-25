import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_velovi import *

method = 'velovi'
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

split1 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split1_v2.h5ad")
vae_split1 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_split1_v2.pt', split1)

split2 = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split2_v2.h5ad")
vae_split2 = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_split2_v2.pt', split2)

total = sc.read_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_total_v2.h5ad")
vae_total = VELOVI.load(data_folder+"v2_"+dataset_long+"/"+method+'/vae_ery_velovi_total_v2.pt', total)

raw = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

#######################################
## add velovi outputs to adata
print_message_with_time("############## Add velovi outputs to adata")

add_velovi_outputs_to_adata(split1, vae_split1)
add_velovi_outputs_to_adata(split2, vae_split2)
add_velovi_outputs_to_adata(total, vae_total)

######################################################
## compute umap
print_message_with_time("############## Compute umap")

compute_umap_ery(split1)
compute_umap_ery(split2)
compute_umap_ery(total)

split1.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split1_outputAdded_v2.h5ad")
split2.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_split2_outputAdded_v2.h5ad")
total.write_h5ad(data_folder+"v2_"+dataset_long+"/"+method+"/adata_"+dataset_short+"_"+method+"_total_outputAdded_v2.h5ad")

#######################################
## plot velocity
print_message_with_time("############## Plot velocity")

plot_velocity_velovi(adata=split1,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,fig_name="split1")
plot_velocity_velovi(adata=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,fig_name="split2")
plot_velocity_velovi(adata=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,fig_name="total")

plot_velocity_velovi(adata=split1,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,fig_name="recompF_split1",recompute=False)
plot_velocity_velovi(adata=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,fig_name="recompF_split2",recompute=False)
plot_velocity_velovi(adata=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder,fig_name="recompF_total",recompute=False)


#######################################
## plot cosine similarity
#compute_cosine_similarity(adata_split1=split1,adata_split2=split2,method='velovi')

#scv.tl.velocity_graph(split1)
#scv.tl.velocity_graph(split2)
#scv.tl.velocity_graph(total
print_message_with_time("############## Plot cosine similarity")

plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_withRef(adata_split1=split1,adata_split2=split2,adata_total=total,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_hist_by_celltype(split1,split2,total,dataset=dataset_short,method=method,fig_folder=fig_folder)

#######################################
## plot velo_conf
print_message_with_time("############## Plot velo_conf")

plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,adata_raw=raw,dataset=dataset_short,method=method,fig_folder=fig_folder)
plot_veloConf_hist(adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,text_x=.7)

#######################################
## plot ptime
# ???

#######################################
## intrinsic uncertainty
"""
save_path = os.path.join(save_dir, "Writeup6_eprs_mg_scvi_ePRS.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    adata2,
    color=["ePRS"],
    frameon=False,
    title="By ePRS status",
    size=5,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')
"""


print_message_with_time("############## Plot intrinsic uncertainty for split1")
compute_intrinisic_uncertainty(adata_in=split1,vae=vae_split1,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split1")

print_message_with_time("############## Plot intrinsic uncertainty for split2")
compute_intrinisic_uncertainty(adata_in=split2,vae=vae_split2,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split2")

print_message_with_time("############## Plot intrinsic uncertainty for total")
compute_intrinisic_uncertainty(adata_in=total,vae=vae_total,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="total")

#######################################
## extrinsic uncertainty
print_message_with_time("############## Plot extrinsic uncertainty for split1")
compute_extrinisic_uncertainty(adata_in=split1,vae=vae_split1,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split1")

print_message_with_time("############## Plot extrinsic uncertainty for split2")
compute_extrinisic_uncertainty(adata_in=split2,vae=vae_split2,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="split2")

print_message_with_time("############## Plot extrinsic uncertainty for total")
compute_extrinisic_uncertainty(adata_in=total,vae=vae_total,adata_raw=raw,dataset=dataset_short,fig_folder=fig_folder,fig_name="total")

#######################################
## permutation score
print_message_with_time("############## Plot permutation score for split1")
compute_permutation_score(adata=split1,vae=vae_split1,dataset=dataset_short,fig_folder=fig_folder,fig_name="split1")

print_message_with_time("############## Plot permutation score for split2")
compute_permutation_score(adata=split2,vae=vae_split2,dataset=dataset_short,fig_folder=fig_folder,fig_name="split2")

print_message_with_time("############## Plot permutation score for split3")
compute_permutation_score(adata=total,vae=vae_total,dataset=dataset_short,fig_folder=fig_folder,fig_name="total")

print_message_with_time("############## All done")
