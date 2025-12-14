import scanpy as sc
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import compute_gene_correlation_between_splits

dataset_long = 'erythroid'
dataset_short = 'ery'
split_seed = 317
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

adata1 = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes_shrunk_disp.h5ad')
adata2 = sc.read(data_folder+'v4_'+dataset_long+'/shrunk_disp/seed'+str(split_seed)+'_'+dataset_short+'_split2_allgenes_shrunk_disp.h5ad')

cor_spliced = compute_gene_correlation_between_splits(adata1.layers['spliced'], adata2.layers['spliced'])
cor_unspliced = compute_gene_correlation_between_splits(adata1.layers['unspliced'], adata2.layers['unspliced'])


import matplotlib.pyplot as plt

plt.clf()
fig, axs = plt.subplots(1, 2, figsize=(12, 4))

# Spliced
axs[0].hist(cor_spliced, bins=30)
axs[0].axvline(np.mean(cor_spliced), color="red", linestyle="--", linewidth=2,
               label=f"Mean = {np.mean(cor_spliced):.3f}")
axs[0].set_xlabel("Correlation")
axs[0].set_ylabel("Frequency")
axs[0].set_title("Spliced gene correlation")
axs[0].legend()

# Unspliced
axs[1].hist(cor_unspliced, bins=30)
axs[1].axvline(np.mean(cor_unspliced), color="red", linestyle="--", linewidth=2,
               label=f"Mean = {np.mean(cor_unspliced):.3f}")
axs[1].set_xlabel("Correlation")
axs[1].set_ylabel("Frequency")
axs[1].set_title("Unspliced gene correlation")
axs[1].legend()

plt.tight_layout()
plt.show()

fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/v4_erythroid/shrunk_disp/seed317/'
plt.savefig(fig_folder+"gene_correlation_shrunk_disp.png", dpi=300, bbox_inches="tight")



