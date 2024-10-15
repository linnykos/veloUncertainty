import anndata as ad

adata_larry_mono = ad.read_h5ad("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry_mono.h5ad")
adata_larry_mono = adata_larry_mono[adata_larry_mono.obs.state_info != "Centroid", :]

adata_larry_neu = ad.read_h5ad("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry_neu.h5ad")
adata_larry_neu = adata_larry_neu[adata_larry_neu.obs.state_info != "Centroid", :]

adata = adata_larry_mono.concatenate(adata_larry_neu)

adata
adata.obs['state_info'].value_counts()
adata.obs['time_info'].value_counts()

adata.write("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry_multilineage.h5ad")