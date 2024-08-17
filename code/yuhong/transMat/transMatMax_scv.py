import scanpy as sc
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_transMat import *

method = 'scv'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

#################################################
#################################################
## erythroid
dataset_long = 'erythroid'
dataset_short = 'ery'
celltype_label = 'celltype'

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

"""
total = sc.read_h5ad(data_folder+"Gastrulation/erythroid_lineage.h5ad")
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/seed317_split1_allgenes.h5ad') # 9815 × 53801
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/seed317_split2_allgenes.h5ad')
gene_names = total.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

S_mat_split1 = split1.layers['spliced'].copy()
U_mat_split1 = split1.layers['unspliced'].copy()
S_mat_split2 = split2.layers['spliced'].copy()
U_mat_split2 = split2.layers['unspliced'].copy()
S_mat_total = total.layers['spliced'].copy()
U_mat_total = total.layers['unspliced'].copy()

def scv_compute_velocity(adata):
    import bbknn
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30) # n_neighbors=30 by default
    ### batch correction
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("Batch correction done!")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 by default
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

scv_compute_velocity(total)
positions_total = [positions_dict[gene] for gene in total.var.index]
total.layers['spliced_original'] = S_mat_total[:,positions_total]
total.layers['unspliced_original'] = U_mat_total[:,positions_total]

scv_compute_velocity(split1)
positions_split1 = [positions_dict[gene] for gene in split1.var.index]
split1.layers['spliced_original'] = S_mat_split1[:,positions_split1] 
split1.layers['unspliced_original'] = U_mat_split1[:,positions_split1]

scv_compute_velocity(split2)
positions_split2 = [positions_dict[gene] for gene in split2.var.index]
split2.layers['spliced_original'] = S_mat_split2[:,positions_split2]
split2.layers['unspliced_original'] = U_mat_split2[:,positions_split2]

total.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')
"""
total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=True,correct_c_same_neighbor=True)
"""
Results:
Number of same celltype prediction = 6693
Proportion of same celltype prediction = 0.7135
Blood progenitors 1: 0.8457
Blood progenitors 2: 0.6309
Erythroid1: 0.7522
Erythroid2: 0.3752
Erythroid3: 0.8752
"""
calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=False,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 6693
Proportion of same celltype prediction = 0.7135
0.8457
0.6309
0.7522
0.3752
0.8752
"""

#################################################
## pancreas
dataset_long = 'pancreas'
dataset_short = 'pan'
celltype_label = 'clusters'

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"
"""
total = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
adata_split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/seed317_split1_allgenes.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/seed317_split2_allgenes.h5ad')
gene_names = total.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()
S_mat_total = total.layers['spliced'].copy()
U_mat_total = total.layers['unspliced'].copy()

## run model
def scv_compute_velocity_pancreas(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

scv_compute_velocity_pancreas(total)
positions_total = [positions_dict[gene] for gene in total.var.index]
total.layers['spliced_original'] = S_mat_total[:,positions_total]
total.layers['unspliced_original'] = U_mat_total[:,positions_total]

scv_compute_velocity_pancreas(adata_split1)
positions_split1 = [positions_dict[gene] for gene in adata_split1.var.index]
adata_split1.layers['spliced_original'] = S_mat_split1[:,positions_split1] 
adata_split1.layers['unspliced_original'] = U_mat_split1[:,positions_split1]

scv_compute_velocity_pancreas(adata_split2)
positions_split2 = [positions_dict[gene] for gene in adata_split2.var.index]
adata_split2.layers['spliced_original'] = S_mat_split2[:,positions_split2]
adata_split2.layers['unspliced_original'] = U_mat_split2[:,positions_split2]

# write data
total.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad')
adata_split1.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
adata_split2.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')
"""

total=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad')
adata_split1=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
adata_split2=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

calculate_transMat_alignment(split1=adata_split1,split2=adata_split2,celltype_label=celltype_label,print_celltype=True,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 2435
Proportion of same celltype prediction = 0.7365
Ductal: 0.9459
Ngn3 low EP: 0.771
Ngn3 high EP: 0.7206
Pre-endocrine: 0.5423
Beta: 0.9433
Alpha: 0.5288
Delta: 0.6377
Epsilon: 0.2606
"""
calculate_transMat_alignment(split1=adata_split1,split2=adata_split2,celltype_label=celltype_label,print_celltype=False,correct_c_same_neighbor=True)
"""
Number of same celltype prediction = 2435
Proportion of same celltype prediction = 0.7365
0.9459
0.771
0.7206
0.5423
0.9433
0.5288
0.6377
0.2606
"""

##################################################
##################################################
dataset_long = "pancreasINC"
dataset_short = "panINC"
celltype_label = 'clusters'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_pancreasINC/'+method+'/'

# pancreasINC_split1_allgenes.h5ad  pancreasINC_split2_allgenes.h5ad  pancreasINC_total_allgenes.h5ad
"""
#raw = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")
total = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_total_allgenes.h5ad')
adata_split1 = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_split1_allgenes.h5ad')
adata_split2 = sc.read_h5ad(data_folder+'v2_pancreasINC/pancreasINC_split2_allgenes.h5ad')
gene_names = total.var.index.copy()
positions_dict = {gene: pos for pos, gene in enumerate(gene_names)}

S_mat_split1 = adata_split1.layers['spliced'].copy()
U_mat_split1 = adata_split1.layers['unspliced'].copy()
S_mat_split2 = adata_split2.layers['spliced'].copy()
U_mat_split2 = adata_split2.layers['unspliced'].copy()
S_mat_total = total.layers['spliced'].copy()
U_mat_total = total.layers['unspliced'].copy()

## run model
def scv_compute_velocity_pancreas(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # previously n_neighbors=10
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

scv_compute_velocity_pancreas(total) 
positions_total = [positions_dict[gene] for gene in total.var.index]
total.layers['spliced_original'] = S_mat_total[:,positions_total]
total.layers['unspliced_original'] = U_mat_total[:,positions_total]

scv_compute_velocity_pancreas(adata_split1)
positions_split1 = [positions_dict[gene] for gene in adata_split1.var.index]
adata_split1.layers['spliced_original'] = S_mat_split1[:,positions_split1] 
adata_split1.layers['unspliced_original'] = U_mat_split1[:,positions_split1]

scv_compute_velocity_pancreas(adata_split2)
positions_split2 = [positions_dict[gene] for gene in adata_split2.var.index]
adata_split2.layers['spliced_original'] = S_mat_split2[:,positions_split2]
adata_split2.layers['unspliced_original'] = U_mat_split2[:,positions_split2]

# write data
total.write_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_scv_total_v2.h5ad')
adata_split1.write_h5ad(data_folder+'v2_'+dataset_long+'/scv/adata_'+dataset_short+'_scv_split1_v2.h5ad')
adata_split2.write_h5ad(data_folder+'v2_'+dataset_long+'/scv/adata_'+dataset_short+'_scv_split2_v2.h5ad')
"""
total=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_scv_total_v2.h5ad')
adata_split1=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/scv/adata_'+dataset_short+'_scv_split1_v2.h5ad')
adata_split2=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/scv/adata_'+dataset_short+'_scv_split2_v2.h5ad')

calculate_transMat_alignment(split1=adata_split1,split2=adata_split2,celltype_label=celltype_label,print_celltype=True,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 2050
Proportion of same celltype prediction = 0.831
Ductal: 0.9513
Ngn3 low EP: 0.8588
Ngn3 high EP: 0.9502
Beta: 0.9707
Alpha: 0.6004
Delta: 0.7
Epsilon: 0.2296
"""

"""
### Results:
Number of same celltype prediction = 2050
Proportion of same celltype prediction = 0.831
0.9513
0.8588
0.9502
0.9707
0.6004
0.7
0.2296
"""
