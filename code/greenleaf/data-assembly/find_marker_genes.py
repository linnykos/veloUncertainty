import anndata as ad
import pandas as pd
import scanpy as sc
from pathlib import Path                
from scipy.io import mmwrite 
import mygene

# from https://github.com/linnykos/veloUncertainty/blob/main/code/greenleaf/create_data/GPC_1data.py
data_folder = '/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
adata = ad.read_h5ad(data_folder+'v4_greenleaf/glf_total_allgenes.h5ad') 
adata.X = adata.layers['spliced'].copy() # we'll do marker genes based on spliced

sc.tl.rank_genes_groups(
    adata,
    groupby="cluster_name",      # use your cluster labels
    method="wilcoxon",           # common choice for scRNA-seq
    n_genes=25,                  # number of genes to rank per cluster
    use_raw=False,               # use adata.X as-is
    key_added="rank_cluster_markers"
)

# 2) Collect top 25 markers per cluster into a dict
top_markers = {}
for cl in adata.obs["cluster_name"].unique():
    df_cl = sc.get.rank_genes_groups_df(
        adata,
        group=cl,
        key="rank_cluster_markers"
    )
    # df_cl is already sorted by the test 'scores' (descending)
    top_markers[cl] = df_cl["names"].head(25).tolist()

# 3) Turn into a tidy DataFrame and/or save to CSV
marker_df = sc.get.rank_genes_groups_df(
    adata,
    group=None,
    key="rank_cluster_markers"
)
# marker_df has columns: 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'group'

#####

# from https://github.com/linnykos/veloUncertainty/blob/main/code/greenleaf/gsea/greenleaf_scvelo_chung_neu4.py
# Initialize
mg = mygene.MyGeneInfo()
# List of all gene IDs from your dataset
ensembl_ids = marker_df['names'].tolist()  # Get all gene IDs from your AnnData object
# Query the mapping
out = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
# Build a mapping: Ensembl ID -> gene symbol
id_to_symbol = {entry['query']: entry.get('symbol', None) for entry in out}
# Print a few mappings to check
for k, v in list(id_to_symbol.items())[:10]:
    print(f"{k} → {v}")

# Map Ensembl IDs to gene symbols
marker_df["gene_symbol"] = [id_to_symbol.get(gene, None) for gene in ensembl_ids]


marker_df.to_csv("/home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/csv/kevin/greenleaf/marker_genes_top25.csv", index=False)

######

n_unique_genes = marker_df["gene_symbol"].nunique()
print(n_unique_genes)

# Get unique, non-missing gene symbols
unique_genes = (
    marker_df["gene_symbol"]
    .dropna()
    .drop_duplicates()   # preserves first-seen order
)

# Write to file: one gene per line, no index, no header, no quotes
unique_genes.to_csv(
    "/home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/csv/kevin/greenleaf/unique_marker_genes.csv",
    index=False,
    header=False
)