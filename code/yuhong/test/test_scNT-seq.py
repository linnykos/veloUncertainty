import dynamo as dyn    

neuron_labeling = dyn.sample_data.scNT_seq_neuron_labeling()
neuron_splicing = dyn.sample_data.scNT_seq_neuron_splicing()  

adata_splicing = neuron_splicing[neuron_labeling.obs_names,[i in neuron_labeling.var_names.values for i in neuron_splicing.var_names.values]]
adata_labeling = neuron_labeling[adata_splicing.obs_names,[i in adata_splicing.var_names.values for i in neuron_labeling.var_names.values]]
adata_splicing.obs = adata_labeling.obs.copy()
"""
>>> adata_splicing
AnnData object with n_obs x n_vars = 3060 x 23918
    obs: 'cell_type', '4sU', 'treatment', 'time', 'cell_name'
    var: 'gene_short_name'
    layers: 'spliced', 'unspliced'
>>> adata_labeling
View of AnnData object with n_obs x n_vars = 3060 x 23918
    obs: 'cell_type', '4sU', 'treatment', 'time', 'cell_name'
    var: 'gene_short_name', 'activity_genes'
    layers: 'new', 'total'
"""

adata_splicing.write_h5ad('data/adata_splicing.h5ad')
adata_labeling.write_h5ad('data/adata_labeling.h5ad')