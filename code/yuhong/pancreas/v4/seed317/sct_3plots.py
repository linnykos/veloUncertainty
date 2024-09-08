split_seed = 317
method = 'sct'
dataset_long = 'pancreas'
dataset_short = 'pan'
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
timestep = df.iloc[np.where(df['pseudotime_corr']==np.max(df['pseudotime_corr']))[0][0]]['time'] # 0.01

adata_prefix = 'adata_'+dataset_short+'_'+method
tnode_prefix = 'tnode_'+dataset_short+'_'+method

"""
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=False)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=False)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=False)
tnode_total = sct.predict.load_model(data_folder+tnode_prefix+'_total_v4.pth')
tnode_split1 = sct.predict.load_model(data_folder+tnode_prefix+'_split1_v4.pth')
tnode_split2 = sct.predict.load_model(data_folder+tnode_prefix+'_split2_v4.pth')

compute_umap(split1, dataset_short)
compute_umap(split2, dataset_short)
compute_umap(total, dataset_short)

sct_seed=615
torch.manual_seed(sct_seed)
random.seed(sct_seed)
np.random.seed(sct_seed)
total.layers['velocity'] = compute_sctour_velocity(tnode_total, timestep=timestep) #+
split1.layers['velocity'] = compute_sctour_velocity(tnode_split1, timestep=timestep) #+
split2.layers['velocity'] = compute_sctour_velocity(tnode_split2, timestep=timestep) #+

total.write(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1.write(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2.write(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

"""
total = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='total',allgenes=False,outputAdded=True)
split1 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split1',allgenes=False,outputAdded=True)
split2 = read_data_v4(dataset_long,dataset_short,method,split_seed,data_version='split2',allgenes=False,outputAdded=True)

############################################
# vector field

plot_vf_umap(adata_in=split1, data_version="split1",data=dataset_short,method=method,fig_folder=fig_folder)
plot_vf_umap(adata_in=split2, data_version="split2",data=dataset_short,method=method,fig_folder=fig_folder)
plot_vf_umap(adata_in=total, data_version="total",data=dataset_short,method=method,fig_folder=fig_folder)

ptime_sct_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)

############################################
## velocity
plot_sct_velocity(adata_in=total,data_version='total',dataset=dataset_short,fig_folder=fig_folder,recompute=True,method='sct',celltype_label=None)
plot_sct_velocity(adata_in=split1,data_version='split1',dataset=dataset_short,fig_folder=fig_folder,recompute=True,method='sct',celltype_label=None)
plot_sct_velocity(adata_in=split2,data_version='split2',dataset=dataset_short,fig_folder=fig_folder,recompute=True,method='sct',celltype_label=None)


############################################
## cosine similarity
plot_cosine_similarity(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
plot_cosine_similarity_hist_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)
plot_cosine_similarity_boxplot_by_celltype(adata_split1=split1,adata_split2=split2,adata_total=total,dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed)

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method) # 1486
c2,n2 = compute_cosine_similarity_union(split1,split2,method) # 2514
np.quantile(c1,[0.,.25,.5,.75,1.]) # [-0.19761002,  0.83260471,  0.88979718,  0.92645803,  0.97551578]
np.quantile(c2,[0.,.25,.5,.75,1.]) # [-0.19190087,  0.81103752,  0.87097793,  0.90436934,  0.95972185]

######################################################
## plot velo_conf
if (not 'velocity_confidence' in total.obs.columns):
    scv.tl.velocity_confidence(total)
    scv.tl.velocity_confidence(split1)
    scv.tl.velocity_confidence(split2)

plot_veloConf_and_cosSim(total,split1,split2,dataset_short,method,fig_folder,split_seed)
plot_veloConf_hist(total,dataset_short,method,fig_folder,split_seed)
plot_velo_conf_boxplot_by_celltype(total,dataset_short,method,fig_folder,split_seed,celltype_label=None)

np.corrcoef(split1.obs['velocity_confidence'],split2.obs['velocity_confidence']) # 0.27310003


######################################################
## ptime
scv.tl.velocity_pseudotime(total,use_velocity_graph=False)
scv.tl.velocity_pseudotime(split1,use_velocity_graph=False)
scv.tl.velocity_pseudotime(split2,use_velocity_graph=False)

plot_pseudotime(adata_in=split1,data_version='split1',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=split2,data_version='split2',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')
plot_pseudotime(adata_in=total,data_version='total',dataset=dataset_short,method=method,fig_folder=fig_folder,split_seed=split_seed,ptime_label='velocity_pseudotime')

ptime_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,time_label='velocity_pseudotime',split_seed=split_seed)
# -0.463



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

latent_time_correlation_scatter_spearman(s1=split1,s2=split2,method=method,dataset=dataset_short,name='split1vs2',xlab='split1',ylab='split2',fig_folder=fig_folder,split_seed=split_seed)
# -0.458

print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='velocity_pseudotime')
print_ptime_corr_by_celltype(split1,split2,total,dataset_short,ptime_label='ptime')

"""
np.corrcoef([c2,total.obs['velocity_confidence'],total.obs['ptime'],total.obs['velocity_pseudotime'],total.obs['latent_time']])
array([[ 1.        ,  0.31381099, -0.06183637,  0.03091665,  0.04335196],
       [ 0.31381099,  1.        , -0.2964452 ,  0.05444068,  0.37584017],
       [-0.06183637, -0.2964452 ,  1.        ,  0.64635625, -0.8548869 ],
       [ 0.03091665,  0.05444068,  0.64635625,  1.        , -0.32460383],
       [ 0.04335196,  0.37584017, -0.8548869 , -0.32460383,  1.        ]])
"""

###################################
# method-selected gene corr
plot_method_gene_corr(split1, split2, method, dataset_short, fig_folder, split_seed)


####################################
# shuffled cosine similarity
v2s_mean,v2s_median = compute_cosine_similarity_shuffled(split1,split2,method=method,seed=1508)
np.round(np.mean(v2s_mean),5) # 
np.round(np.mean(v2s_median),5) # 

np.var(v2s_mean) # 2.518031761592001e-05
np.var(v2s_median) # 6.200052935986668e-06
# paired: 0.8446 0.871
# shuffled: 0.11892 -0.02414
### shuffled mean and median have variance approximately 0 among 100 iterations