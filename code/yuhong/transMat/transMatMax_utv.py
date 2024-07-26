import scanpy as sc
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import *
from v2_functions_transMat import *
from v2_functions_transMat import calculate_transMat_alignment

method = 'utv'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

#################################################
#################################################
## erythroid
dataset_long = 'erythroid'
dataset_short = 'ery'
celltype_label = 'celltype'

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

def utv_compute_umap(adata):
    import scvelo as scv
    import bbknn
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("************ batch correction done ************")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # used to be 10
    sc.tl.umap(adata)

utv_compute_umap(total)
utv_compute_umap(split1)
utv_compute_umap(split2)

calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=False,correct_c_same_neighbor=True)
calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=True,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 7136
Proportion of same celltype prediction = 0.7607
Blood progenitors 1: 0.8955
Blood progenitors 2: 0.6575
Erythroid1: 0.5624
Erythroid2: 0.9277
Erythroid3: 1.0
"""

#################################################
#################################################
## pancreas
dataset_long = 'pancreas'
dataset_short = 'pan'
celltype_label = 'clusters'

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"

total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad')

sc.tl.umap(total)
sc.tl.umap(split1)
sc.tl.umap(split2)

calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=True,correct_c_same_neighbor=False,use_negative_cosines=True)
calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=False,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 2214
Proportion of same celltype prediction = 0.6329
Ductal: 0.4729
Ngn3 low EP: 0.5709
Ngn3 high EP: 0.749
Pre-endocrine: 0.4778
Beta: 0.7642
Alpha: 0.9396
Delta: 0.1143
Epsilon: 0.7042
"""
c_tmp = 0
for i in np.where(list(split1.obs['clusters']=='Delta'))[0]: 
    pred_idx1 = np.where(t1[i].todense()==np.max(t1[i]))[1]
    celltype1 = split1.obs['clusters'][pred_idx1].to_numpy()[0]
    pred_idx2 = np.where(t2[i].todense()==np.max(t2[i]))[1]
    celltype2 = split2.obs['clusters'][pred_idx2].to_numpy()[0]
    if celltype1==celltype2: c_tmp+=1
    print(celltype1, celltype2)
# many split1 Delta cells point to themselves, i.e. not confident in predicting

"""
>>> split1.obsp['connectivities']
<3696x3696 sparse matrix of type '<class 'numpy.float32'>'
        with 151814 stored elements in Compressed Sparse Row format>
>>> split2.obsp['connectivities']
<3696x3696 sparse matrix of type '<class 'numpy.float32'>'
        with 151714 stored elements in Compressed Sparse Row format>
>>> np.sum(split1.obsp['connectivities'][0].todense()>0)
81
>>> np.sum(split1.obsp['connectivities'][1].todense()>0)
38
>>> np.sum(split1.obsp['connectivities'][2].todense()>0)
30
>>> np.sum(split1.obsp['connectivities'][3].todense()>0)
38
>>> np.sum(split1.obsp['connectivities'][4].todense()>0)
29
"""
def compute_umap(adata):
    import scvelo as scv
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # used to be n_neighbors=10
    sc.tl.umap(adata)
    scv.tl.velocity_graph(adata)

compute_umap(split1)
compute_umap(split2)
"""
>>> split1.obsp['connectivities']
<3696x3696 sparse matrix of type '<class 'numpy.float32'>'
        with 76602 stored elements in Compressed Sparse Row format>
>>> np.sum(split1.obsp['connectivities'][0].todense()>0)
33
>>> np.sum(split1.obsp['connectivities'][1].todense()>0)
18
>>> np.sum(split1.obsp['connectivities'][2].todense()>0)
14
>>> np.sum(split1.obsp['connectivities'][3].todense()>0)
20
>>> np.sum(split1.obsp['connectivities'][4].todense()>0)
14
"""

calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=True,correct_c_same_neighbor=True,use_negative_cosines=False)
calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=False,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 2071
Proportion of same celltype prediction = 0.6264
Ductal: 0.611
Ngn3 low EP: 0.5649
Ngn3 high EP: 0.7428
Pre-endocrine: 0.4965
Beta: 0.563
Alpha: 0.8145
Delta: 0.3913
Epsilon: 0.6901
"""

#################################################
#################################################
## pancreasINC
import numpy as np
import anndata as ad
dataset_long = 'pancreasINC'
dataset_short = 'panINC'
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"
celltype_label = 'clusters'

#total = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v2.h5ad') # 
split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 
#raw = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")

compute_umap(split1)
compute_umap(split2)

calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=True,correct_c_same_neighbor=True)
calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=False,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 1681
Proportion of same celltype prediction = 0.6814
Ductal: 0.6059
Ngn3 low EP: 0.5534
Ngn3 high EP: 0.9502
Beta: 0.6117
Alpha: 0.8774
Delta: 0.6429
Epsilon: 0.5852
"""

#################################################
#################################################
## LARRY
dataset_long = 'larry'
dataset_short = 'larry'
fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_"+dataset_long+"/"+method+"/"
celltype_label = 'state_info'

split1 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v2.h5ad') # 49302 x 2000
split2 = sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v2.h5ad') # 49302 x 2000

import matplotlib as mpl
def rgb2hex(rgb):
    r = int(rgb[0]*255)
    g = int(rgb[1]*255)
    b = int(rgb[2]*255)
    return "#{:02x}{:02x}{:02x}".format(r,g,b)

split1.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]
split2.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[40:315:25]]

compute_umap(split1)
compute_umap(split2)

calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=True,correct_c_same_neighbor=True)
calculate_transMat_alignment(split1,split2,celltype_label,print_celltype=False,correct_c_same_neighbor=True)
"""
### Results:
Number of same celltype prediction = 33069
Proportion of same celltype prediction = 0.8344
Baso: 0.9026
Ccr7_DC: 0.875
Eos: 0.5833
Erythroid: 0.8219
Lymphoid: 0.8812
Mast: 0.9517
Meg: 0.9785
Monocyte: 0.5865
Neutrophil: 0.6465
Undifferentiated: 0.9701
pDC: 0.6122
"""