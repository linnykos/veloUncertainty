import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_larry/" 

# correlation between splits
split1_allgenes = sc.read_h5ad(data_folder+'v2_larry/larry_split1_allgenes.h5ad')
split2_allgenes = sc.read_h5ad(data_folder+'v2_larry/larry_split2_allgenes.h5ad')
split1_allgenes.layers['spliced_original'] = split1_allgenes.layers['spliced']
split1_allgenes.layers['unspliced_original'] = split1_allgenes.layers['unspliced']
split2_allgenes.layers['spliced_original'] = split2_allgenes.layers['spliced']
split2_allgenes.layers['unspliced_original'] = split2_allgenes.layers['unspliced']

plot_gene_correlation_between_splits(adata1=split1_allgenes,adata2=split2_allgenes,fig_path='corr_between_splits_allgenes.png',fig_folder=fig_folder)

# mark housekeeping genes
housekeeping_genes = pd.read_csv("/home/users/y2564li/kzlinlab/data/genelists/housekeeping/HRT_Atlas/Housekeeping_GenesMouse_formatted.csv",header=None)
housekeeping_genes = housekeeping_genes[0] # 3277 genes

hkgenes_common = split1_allgenes.var.index.intersection(housekeeping_genes)

positions_dict = {gene: pos for pos, gene in enumerate(split1_allgenes.var.index)}
positions = [positions_dict[gene] for gene in hkgenes_common]

cor_spliced = compute_gene_correlation_between_splits(split1_allgenes.layers['spliced'],split2_allgenes.layers['spliced'])
cor_unspliced = compute_gene_correlation_between_splits(split1_allgenes.layers['unspliced'],split2_allgenes.layers['unspliced'])
Ngenes_spliced = len(cor_spliced[~np.isnan(cor_spliced)])
Ngenes_unspliced = len(cor_unspliced[~np.isnan(cor_unspliced)])
colors = ['royalblue']*len(cor_spliced)
for i in positions: 
    colors[i] = 'indianred'

colors = np.array(colors)

# plot allgenes, with housekeeping genes marked
## spliced
plt.scatter(range(Ngenes_spliced), cor_spliced[~np.isnan(cor_spliced)],color=colors[~np.isnan(cor_spliced)],alpha=0.4)
plt.title("Correlation of gene expr between splits (spliced, housekeeping genes marked), N="+str(Ngenes_spliced))
plt.ylim(-.1, 1)
plt.savefig(fig_folder+'corr_between_splits_allgenes_markHKgenes_spliced.png') 
plt.clf()
## unspliced
plt.scatter(range(Ngenes_unspliced), cor_unspliced[~np.isnan(cor_unspliced)],color=colors[~np.isnan(cor_unspliced)],alpha=0.4)
plt.title("Correlation of gene expr between splits (unspliced, housekeeping genes marked), N="+str(Ngenes_unspliced))
plt.ylim(-.1, 1)
plt.savefig(fig_folder+'corr_between_splits_allgenes_markHKgenes_unspliced.png') 
plt.clf()
## spliced and unspliced together
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
axes[0].scatter(range(Ngenes_spliced), cor_spliced[~np.isnan(cor_spliced)],color=colors[~np.isnan(cor_spliced)],alpha=0.4)
axes[0].set_title("Correlation of gene expr between splits (spliced, housekeeping genes marked), N="+str(Ngenes_spliced))
axes[0].set_xlabel("Genes")
axes[0].set_ylabel("Correlation")
axes[1].scatter(range(Ngenes_unspliced), cor_unspliced[~np.isnan(cor_unspliced)],color=colors[~np.isnan(cor_unspliced)],alpha=0.4)
axes[1].set_title("Correlation of gene expr between splits (unspliced, housekeeping genes marked), N="+str(Ngenes_unspliced))
axes[1].set_xlabel("Genes")
axes[1].set_ylabel("Correlation")
plt.tight_layout() # Adjusting layout to avoid overlap
plt.savefig(fig_folder+'corr_between_splits_allgenes_markHKgenes.png') 
plt.clf()
###




