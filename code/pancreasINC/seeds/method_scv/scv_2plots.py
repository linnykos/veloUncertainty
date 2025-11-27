import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

def scv_panINC_plots(split_seed):
       print('############### seed='+str(split_seed)+' starts! ###############')
       dataset_long = 'pancreasINC'
       dataset_short = 'panINC'
       method = 'scv'
       #data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
       fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
       celltype_label = get_celltype_label(dataset_short)
       print(str(split_seed)+': read data')
       total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
       split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
       split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
       # method-selected gene corr
       print(str(split_seed)+': gene correlation plots')
       plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)
       ## velocity
       print(str(split_seed)+': velocity plots')
       plot_velocity_scv_utv(adata_in=split1,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
       plot_velocity_scv_utv(adata_in=split2,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
       ## cosine similarity
       print(str(split_seed)+': cosine similarity')
       plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
       plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
       plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
       plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
       ## confidence
       print(str(split_seed)+': velocity confidence')
       scv.tl.velocity_confidence(total)
       scv.tl.velocity_confidence(split1)
       scv.tl.velocity_confidence(split2)
       plot_veloConf_and_cosSim(adata_total=total,adata_split1=split1,adata_split2=split2,dataset=dataset_short,method=method,fig_folder=fig_folder, split_seed=split_seed)
       plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
       plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)
       print('correlation of confidence between splits ='+str(np.round(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']),5))) # 0.28853776
       # shuffled cosine similarity
       print(str(split_seed)+': shuffled cosine similarity')
       v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
       print('mean of shuffled cosine similarity mean='+str(np.round(np.mean(v2s_mean),5)) )
       print('mean of shuffled cosine similarity median ='+str(np.round(np.mean(v2s_median),5)) )
       # cosine similarity, intersection and union of genes
       print(str(split_seed)+':cosine similarity, intersected genes')
       c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
       print(np.quantile(c1,[0.,.25,.5,.75,1.]) )
       np.round(np.median(c1),4) 
       np.round(np.mean(c1),4) 
       np.round(np.var(c1),4)
       print(str(split_seed)+':cosine similarity, unioned genes')
       c2,n2 = compute_cosine_similarity_union(split1,split2,method)
       np.quantile(c2,[0.,.25,.5,.75,1.]) 
       np.round(np.median(c2),4) 
       np.round(np.mean(c2),4) 
       np.round(np.var(c2),4) 
       print('############### seed='+str(split_seed)+' all done! ###############')

scv_panINC_plots(320)
scv_panINC_plots(323)
scv_panINC_plots(326)
scv_panINC_plots(329)
