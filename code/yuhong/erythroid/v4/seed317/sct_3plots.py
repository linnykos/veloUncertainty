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

df = pd.read_csv(data_folder+dataset_short+'_sct_velo_timestep.csv')
timestep = df.iloc[np.where(df['pseudotime_corr']==np.max(df['pseudotime_corr']))[0][0]]['time'] # 0.19

adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method

"""
total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4.h5ad') # 
tnode_total = sct.predict.load_model(data_folder+tnode_prefix+'_total_v4.pth')
tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')

sct_add_state_info_colors(split1)
sct_add_state_info_colors(split2)
sct_add_state_info_colors(total)

"""
total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 


############################################
# vector field

plot_vf_umap(adata_in=split1, data_version="split1",data=dataset_short,method=method,fig_folder=fig_folder)
plot_vf_umap(adata_in=split2, data_version="split2",data=dataset_short,method=method,fig_folder=fig_folder)
plot_vf_umap(adata_in=total, data_version="total",data=dataset_short,method=method,fig_folder=fig_folder)

ptime_sct_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)

############################################
# velocity
"""
total.layers['velocity'] = compute_sctour_velocity(tnode_total, timestep=timestep) #+
split1.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=timestep) #-
split2.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=timestep) #+
print('Velocity computed')

get_umap_sct(total)
get_umap_sct(split1)
get_umap_sct(split2)
print('UMAP computed')

total.write_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1.write_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2.write_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 
"""
plot_sct_velocity(adata_in=split1,data_version='split1',dataset=dataset_short,fig_folder=fig_folder)
plot_sct_velocity(adata_in=split2,data_version='split2',dataset=dataset_short,fig_folder=fig_folder)
plot_sct_velocity(adata_in=total,data_version='total',dataset=dataset_short,fig_folder=fig_folder)

######################################################
## plot cosine similarity
#c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
#c2,n2 = compute_cosine_similarity_union(split1,split2,method)
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

"""
>>> np.corrcoef([c2,total.obs['ptime'],total.obs['velocity_confidence']])
array([[ 1.        ,  0.39571204,  0.39549537],
       [ 0.39571204,  1.        , -0.2266979 ],
       [ 0.39549537, -0.2266979 ,  1.        ]])
>>> np.corrcoef([c2,total.obs['velocity_pseudotime'],total.obs['velocity_confidence']])
array([[ 1.        ,  0.37199716,  0.39549537],
       [ 0.37199716,  1.        , -0.17789574],
       [ 0.39549537, -0.17789574,  1.        ]])
"""
######################################################
## plot velo_conf
scv.tl.velocity_confidence(total)
scv.tl.velocity_confidence(split1)
scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed)

print(np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence'])) # 0.13712638

######################################################
## ptime
scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)

plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)


plot_pseudotime_diffusion(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime_diffusion(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime_diffusion(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')

ptime_diffusion_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)


print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='velocity_pseudotime')
print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='ptime')

"""
>>> print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='velocity_pseudotime')
Blood progenitors 1 0.9098
Blood progenitors 2 0.8911
Erythroid1 0.6146
Erythroid2 0.4917
Erythroid3 0.5882
>>> print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='ptime')
Blood progenitors 1 0.9308
Blood progenitors 2 0.9325
Erythroid1 0.7669
Erythroid2 0.5862
Erythroid3 0.6855
"""

if not 'latent_time' in split1.obs.columns:
    scv.tl.recover_dynamics(total,n_jobs=8)
    scv.tl.latent_time(total)
    plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    scv.tl.recover_dynamics(split1,n_jobs=8)
    scv.tl.latent_time(split1)
    plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
    scv.tl.recover_dynamics(split2,n_jobs=8)
    scv.tl.latent_time(split2)
    plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')


#plot_latent_time(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
#plot_latent_time(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')
#plot_latent_time(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,celltype_label='state_info')

###################################
# method-selected gene corr
plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)


####################################
# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
np.round(np.mean(v2s_mean),5) # 0.47158
np.round(np.mean(v2s_median),5) # 0.47913

# paired: 0.6231 0.6483
# shuffled: 0.47158 0.47913
### shuffled mean and median have variance approximately 0 among 100 iterations

####################################
cos_sim_vf = np.diag(cosine_similarity(split1.obsm['X_VF'],split2.obsm['X_VF']))
[np.round(np.median(cos_sim_vf),4), np.round(np.mean(cos_sim_vf),4)]

