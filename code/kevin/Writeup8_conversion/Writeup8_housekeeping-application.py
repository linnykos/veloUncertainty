import anndata as ad
import pandas as pd

adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/v3_erythroid/sct/adata_ery_sct_total.h5ad")
housekeeping_genes = pd.read_csv("/home/users/kzlin/kzlinlab/data/genelists/housekeeping/HRT_Atlas/Housekeeping_GenesMouse_formatted.csv", 
                                 header=None)
housekeeping_genes = housekeeping_genes[0]

housekeeping_genes = adata.var_names.intersection(housekeeping_genes)

print(housekeeping_genes)
print(housekeeping_genes.size)