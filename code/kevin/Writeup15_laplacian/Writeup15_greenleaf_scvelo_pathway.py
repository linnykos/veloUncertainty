# https://github.com/linnykos/veloUncertainty/blob/main/veloUncertainty/v4_functions_transMat.py
# cellrank virtual environment

import scvelo as scv
import numpy as np
import collections

import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

import anndata as ad
import numpy as np
import scipy.sparse as sp
import mygene
import pandas as pd
import gseapy as gp
from biomart import BiomartServer    # or mygene, ensembl REST, etc.
import mygene

# Load the AnnData object
adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/seed317/scv/adata_glf_scv_total_v4.h5ad")

v = np.abs(adata.layers['velocity']).mean(axis=0)         # s_g
gene_rank = pd.Series(v, index=adata.var_names).sort_values(ascending=False)

mg = mygene.MyGeneInfo()
mapping = mg.querymany(gene_rank.index.tolist(), scopes="ensembl.gene",
                       fields="symbol", species="human")
ens2sym = {d["query"]: d.get("symbol") for d in mapping if "symbol" in d}

# replace the index with symbols and drop missing ones
gene_rank_sym = gene_rank.rename(index=ens2sym).dropna()
gene_rank_sym = gene_rank_sym[~gene_rank_sym.index.duplicated()]
gene_rank_sym = gene_rank_sym.sort_values(ascending=False)

res = gp.prerank(
        rnk        = gene_rank_sym,                      # your Series
        gene_sets  = "GO_Biological_Process_2025",   # or another library
        organism   = "Human",                        # tell gseapy the species
        min_size   = 5,
        max_size   = 500,
        outdir     = None,                           # don’t write files
).res2d
res[['Term','NES','FDR q-val']][0:50]

# Assuming your result DataFrame is called `res` and the column is named 'FDR q-val'
significant_pathways = res[res['FDR q-val'] < 0.1]

# View top 10 by NES
print(significant_pathways.sort_values('NES', ascending=False).head(10))