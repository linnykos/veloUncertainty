import scanpy as sc
import scvelo as scv
import sctour as sct
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_sct import *


split_seed = 317
method = 'sct'
dataset_long = 'larry'
dataset_short = 'larry'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"
adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method

"""
raw = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'+"v4_larry/larry.h5ad") # n_obs × n_vars = 49302 × 23420

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4.h5ad') # 
tnode_total = sct.predict.load_model(data_folder+tnode_prefix+'_total_v4.pth')
tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')

total.obsm['X_umapOriginal'] = raw.obsm['X_emb'].copy()
split1.obsm['X_umapOriginal'] = raw.obsm['X_emb'].copy()
split2.obsm['X_umapOriginal'] = raw.obsm['X_emb'].copy()

sct_add_state_info_colors(split1)
sct_add_state_info_colors(split2)
sct_add_state_info_colors(total)
print('Read data done')
"""

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

############################################
# vector field

plot_vf_umap(adata_in=split1, data_version="split1",data=dataset_short,method=method,fig_folder=fig_folder,celltype_label='state_info')
plot_vf_umap(adata_in=split2, data_version="split2",data=dataset_short,method=method,fig_folder=fig_folder,celltype_label='state_info')
plot_vf_umap(adata_in=total, data_version="total",data=dataset_short,method=method,fig_folder=fig_folder,celltype_label='state_info')

ptime_sct_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

############################################
"""
timestep = .01
total.layers['velocity'] = compute_sctour_velocity(tnode_total, timestep=timestep)
split1.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=timestep)
split2.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=timestep)
print('Velocity computed')

get_umap_sct(total)
get_umap_sct(split1)
get_umap_sct(split2)
print('UMAP computed')
"""

"""
total.write_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1.write_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2.write_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
"""

"""
plot_sct_velocity(adata_in=split1,data_version='split1',dataset=dataset_short,fig_folder=fig_folder,celltype_label='state_info')
plot_sct_velocity(adata_in=split2,data_version='split2',dataset=dataset_short,fig_folder=fig_folder,celltype_label='state_info')
plot_sct_velocity(adata_in=total,data_version='total',dataset=dataset_short,fig_folder=fig_folder,celltype_label='state_info')

######################################################
## plot cosine similarity
#c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
#c2,n2 = compute_cosine_similarity_union(split1,split2,method)
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

######################################################
## plot velo_conf
if not 'velocity_confidence' in total.obs.columns:
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,raw,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)

print(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']))
"""

######################################################
## ptime
if not 'velocity_pseudotime' in split1.obs.columns:
    scv.tl.velocity_pseudotime(total)
    scv.tl.velocity_pseudotime(split1)
    scv.tl.velocity_pseudotime(split2)

"""
plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info',ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed,celltype_label='state_info')
"""
if not 'latent_time' in split1.obs.columns:
    scv.tl.recover_dynamics(total,n_jobs=8)
    scv.tl.latent_time(total)
    plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    scv.tl.recover_dynamics(split1,n_jobs=8)
    scv.tl.latent_time(split1)
    plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
    scv.tl.recover_dynamics(split2,n_jobs=8)
    scv.tl.latent_time(split2)
    plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')



#plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
#plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
#plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
