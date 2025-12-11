import scanpy as sc
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *

dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'


split_seed = 317
grid_seed = 227

total_mark = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/glf_total_mark_genes.h5ad')
split1_mark = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split1_mark_genes.h5ad')
split2_mark = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'_'+dataset_short+'_split2_mark_genes.h5ad')


method = 'scv'
print_message_with_time('################## start velocity') 
scv_compute_velocity(total_mark, dataset_short) 
scv_compute_velocity(split1_mark, dataset_short) 
scv_compute_velocity(split2_mark, dataset_short) 
print_message_with_time('################## end velocity') 

total_mark.write_h5ad(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_total_mark_genes.h5ad')
split1_mark.write_h5ad(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_split1_mark_genes.h5ad')
split2_mark.write_h5ad(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_'+method+'_split2_mark_genes.h5ad')


### make plots
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *


method = 'scv_mark'
celltype_label = 'cluster_name'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total_mark.obsm['X_umapOriginal'] = total_mark.obsm['X_umap_greenleaf'].copy()
split1_mark.obsm['X_umapOriginal'] = total_mark.obsm['X_umap_greenleaf'].copy()
split2_mark.obsm['X_umapOriginal'] = total_mark.obsm['X_umap_greenleaf'].copy()

## velocity
plot_velocity_scv_utv(adata_in=total_mark,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
plot_velocity_scv_utv(adata_in=split1_mark,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
plot_velocity_scv_utv(adata_in=split2_mark,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)

## cosine similarity
plot_cosine_similarity(adata_split1=split1_mark,adata_split2=split2_mark,adata_total=total_mark,dataset=dataset_short,method=method,
					   fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1_mark, split2_mark, total_mark, dataset_short, method, fig_folder, split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1_mark,adata_split2=split2_mark,adata_total=total_mark,dataset=dataset_short,method=method,
										fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1_mark,adata_split2=split2_mark,adata_total=total_mark,dataset=dataset_short,method=method,
											fig_folder=fig_folder,split_seed=split_seed,celltype_label=celltype_label)

## confidence
print('############### velocity confidence')
scv.tl.velocity_confidence(total_mark)
scv.tl.velocity_confidence(split1_mark)
scv.tl.velocity_confidence(split2_mark)

plot_veloConf_and_cosSim(adata_total=total_mark, adata_split1=split1_mark, adata_split2=split2_mark, dataset=dataset_short, method=method,
						 fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total_mark, dataset_short, method, fig_folder, split_seed)
plot_velo_conf_boxplot_by_celltype(total_mark, dataset_short, method, fig_folder, split_seed, celltype_label=celltype_label)

#######################################
## control
gene_set_name = 'markctrl'
total_ctrl = sc.read(data_folder+'v4_greenleaf/glf_genes_mark/glf_total_'+gene_set_name+str(grid_seed)+'.h5ad')
split1_ctrl = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/glf_split1_'+gene_set_name+str(grid_seed)+'.h5ad')
split2_ctrl = sc.read(data_folder+'v4_'+dataset_long+'/glf_genes_mark/glf_split2_'+gene_set_name+str(grid_seed)+'.h5ad')

print_message_with_time('################## start velocity') 
scv_compute_velocity(total_ctrl, dataset_short) 
scv_compute_velocity(split1_ctrl, dataset_short) 
scv_compute_velocity(split2_ctrl, dataset_short) 
print_message_with_time('################## end velocity') 

total_ctrl.write_h5ad(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_scv_total_markctrl_genes.h5ad')
split1_ctrl.write_h5ad(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_scv_split1_markctrl_genes.h5ad')
split2_ctrl.write_h5ad(data_folder+'v4_'+dataset_long+'/glf_genes_mark/seed'+str(split_seed)+'/scv/adata_'+dataset_short+'_scv_split2_markctrl_genes.h5ad')

### make plots
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

method = 'scv_markctrl'
celltype_label = 'cluster_name'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total_ctrl.obsm['X_umapOriginal'] = total_ctrl.obsm['X_umap_greenleaf'].copy()
split1_ctrl.obsm['X_umapOriginal'] = total_ctrl.obsm['X_umap_greenleaf'].copy()
split2_ctrl.obsm['X_umapOriginal'] = total_ctrl.obsm['X_umap_greenleaf'].copy()

## velocity
plot_velocity_scv_utv(adata_in=total_ctrl,fig_folder=fig_folder,data_version='total',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
plot_velocity_scv_utv(adata_in=split1_ctrl,fig_folder=fig_folder,data_version='split1',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)
plot_velocity_scv_utv(adata_in=split2_ctrl,fig_folder=fig_folder,data_version='split2',dataset=dataset_short,method=method,split_seed=split_seed,celltype_label=celltype_label)

## cosine similarity
plot_cosine_similarity(adata_split1=split1_ctrl, adata_split2=split2_ctrl, adata_total=total_ctrl, dataset=dataset_short, method=method, fig_folder=fig_folder, split_seed=split_seed)
plot_cosine_similarity_withRef(split1_ctrl, split2_ctrl, total_ctrl, dataset_short, method, fig_folder, split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1_ctrl, adata_split2=split2_ctrl, adata_total=total_ctrl, dataset=dataset_short, method=method,
										fig_folder=fig_folder, split_seed=split_seed, celltype_label=celltype_label)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1_ctrl, adata_split2=split2_ctrl, adata_total=total_ctrl, dataset=dataset_short, method=method,
											fig_folder=fig_folder, split_seed=split_seed, celltype_label=celltype_label)

## confidence
print('############### velocity confidence')
scv.tl.velocity_confidence(total_ctrl)
scv.tl.velocity_confidence(split1_ctrl)
scv.tl.velocity_confidence(split2_ctrl)

plot_veloConf_and_cosSim(adata_total=total_ctrl, adata_split1=split1_ctrl, adata_split2=split2_ctrl, dataset=dataset_short, method=method,
						fig_folder=fig_folder, split_seed=split_seed)
plot_veloConf_hist(total_ctrl, dataset_short, method, fig_folder, split_seed)
plot_velo_conf_boxplot_by_celltype(total_ctrl, dataset_short, method, fig_folder, split_seed, celltype_label=celltype_label)



