import anndata as ad
import pandas as pd
import scanpy as sc
from pathlib import Path                
from scipy.io import mmwrite 
import mygene
import seaborn as sns 

# from https://github.com/linnykos/veloUncertainty/blob/main/code/greenleaf/create_data/GPC_1data.py
data_folder = '/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
adata = ad.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 

# from https://github.com/linnykos/veloUncertainty/blob/main/code/greenleaf/gsea/greenleaf_scvelo_chung_neu4.py
# Initialize
mg = mygene.MyGeneInfo()
# List of all gene IDs from your dataset
ensembl_ids = adata.var_names.tolist()  # Get all gene IDs from your AnnData object
# Query the mapping
out = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
# Build a mapping: Ensembl ID -> gene symbol
id_to_symbol = {entry['query']: entry.get('symbol', None) for entry in out}
# Print a few mappings to check
for k, v in list(id_to_symbol.items())[:10]:
    print(f"{k} → {v}")

# Get gene names from adata
gene_names = adata.var_names.tolist()
# Map Ensembl IDs to gene symbols
gene_symbols = [id_to_symbol.get(gene, None) for gene in gene_names]

############

# gene_symbols is your list, same length/order as adata.var_names
# Step 1: make a Series so we can use .duplicated() nicely
gene_symbols_ser = pd.Series(gene_symbols, index=adata.var_names, name="gene_symbol")

# Step 2: build a boolean mask of which genes to KEEP
#  - drop None
#  - drop duplicated symbols (keep the first occurrence)
mask_not_none  = gene_symbols_ser.notna()
mask_unique    = ~gene_symbols_ser.duplicated(keep="first")
mask_keep      = mask_not_none & mask_unique

print(f"Original n_vars: {adata.n_vars}")
print(f"Keeping {mask_keep.sum()} genes after removing None and duplicates")

# Step 3: subset adata to only those genes
adata = adata[:, mask_keep].copy()

# Step 4: store both Ensembl IDs and gene symbols in .var,
#         and make var_names = gene symbols
adata.var["ensembl_id"]  = gene_symbols_ser.index[mask_keep]      # original IDs
adata.var["gene_symbol"] = gene_symbols_ser[mask_keep].values     # unique symbols
adata.var_names = adata.var["gene_symbol"]

############

# from https://github.com/linnykos/ADRC_workshop_2025/blob/main/Day5/seaad-conversion_microglia-pfc_part1.ipynb
# This is the raw count data
adata.X = adata.layers['spliced'].copy() # we'll do marker genes based on spliced counts
adata.X[11:20,21:30].toarray() 

# remove every entry stored in .uns
adata.uns.clear()          # empties the dict in-place
adata

# remove every matrix stored in .layers
adata.layers.clear()      # empties the mapping in-place
adata

# delete all connectivities / distance matrices
adata.obsp.clear()        # empties the mapping in-place
adata


