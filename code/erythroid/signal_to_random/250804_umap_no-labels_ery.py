# scv
import scvelo as scv
import scanpy as sc
import bbknn
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

dataset_long = 'erythroid'
dataset_short = 'ery'
method = 'scv'
split_seed=317
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_scv_total_v4_outputAdded.h5ad')

celltype_label = get_celltype_label(dataset_short)
data_method = dataset_short+"_"+method
data_version = 'total'
total.obsm['X_umap'] = total.obsm['X_umapOriginal'].copy()
scv.pl.velocity_embedding_stream(total, basis='umap',color=celltype_label,recompute=True, legend_loc='none', arrow_size=2,
                                 title=dataset_short+'+'+method, dpi=300, 
                                 save=fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_nolabel.png")    



# utv
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import plot_velocity_scv_utv
from v4_functions import *

split_seed=317
dataset_long = 'erythroid'
dataset_short = 'ery'
method_prefix = 'utv'
method = method_prefix #+ '_' + gene_set_name
celltype_label = get_celltype_label(dataset_short)
data_method = dataset_short+"_"+method
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_ery_utv_total_v4_outputAdded.h5ad')

data_version = 'total'
total.obsm['X_umap'] = total.obsm['X_umapOriginal'].copy()
scv.pl.velocity_embedding_stream(total, basis='umap',color=celltype_label,recompute=True, legend_loc='none', arrow_size=2,
                                 title=dataset_short+'+'+method, dpi=300,
                                 save=fig_folder+"velocity/"+data_method+"_total_"+"umapOriginal_nolabel.png")    

# sct
split_seed = 317
method = 'sct'
dataset_long = 'erythroid'
dataset_short = 'ery'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_sct import *
from v4_functions import *

adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method
total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 

celltype_label = get_celltype_label(dataset_short)
data_method = dataset_short+"_"+method
data_version = 'total'
total.obsm['X_umap'] = total.obsm['X_umapOriginal'].copy()
scv.pl.velocity_embedding_stream(total, basis='umap',color=celltype_label,recompute=True, legend_loc='none', arrow_size=2, 
                                 dpi=300, title=dataset_short+'+'+method,
                                 save=fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_nolabel.png")    

# velovi
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_velovi import *

split_seed = 317
method = 'velovi'
dataset_short = 'ery'
dataset_long = 'erythroid'

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)

celltype_label = get_celltype_label(dataset_short)
data_method = dataset_short+"_"+method
data_version = 'total'
total.obsm['X_umap'] = total.obsm['X_umapOriginal'].copy()
scv.pl.velocity_embedding_stream(total, basis='umap',color=celltype_label,recompute=True,legend_loc='none', arrow_size=2,
                                 dpi=300,title=dataset_short+'+'+method,
                                 save=fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_nolabel.png")    


# without preprocessing
method = 'velovi_woprep'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)

total.obsm['X_umap'] = total.obsm['X_umapOriginal'].copy()
data_method = dataset_short+"_"+method
scv.pl.velocity_embedding_stream(total, basis='umap',color=celltype_label,recompute=True,legend_loc='none', arrow_size=2, 
                                 dpi=300,title=dataset_short+'+'+method,
                                 save=fig_folder+"velocity/"+data_method+"_"+data_version+'_'+"umapOriginal_nolabel.png")    

