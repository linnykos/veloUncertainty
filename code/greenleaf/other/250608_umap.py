dataset_long = 'greenleaf'
dataset_short = 'glf'

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_sct import *

fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'
data = sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/glf_total_allgenes.h5ad')
data.obsm['X_umap'] = data.obsm['X_umap_greenleaf']
scv.settings.figdir = fig_folder
scv.settings.set_figure_params(dpi=150, figsize=(8, 6))

scv.pl.umap(data, color='cluster_name', legend_loc='right margin', save='X_umap_greenleaf_labels_out.png', size=5, alpha=.8)
scv.pl.umap(data, color='cluster_name', legend_loc='right margin')

scv.settings.set_figure_params(dpi=150, figsize=(6,4))
scv.pl.umap(data, color='cluster_name', legend_loc='lower center', save='X_umap_greenleaf_labels_out2.png', size=5, alpha=.8)



import matplotlib.pyplot as plt
scv.pl.umap(data, color='cluster_name', legend_loc='none', save='X_umap_greenleaf_labels_out.png', size=5, alpha=.8)
plt.legend(*plt.gca().get_legend_handles_labels(),
           bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig(fig_folder+'X_umap_greenleaf_labels_out.png')

import matplotlib.pyplot as plt
scv.settings.set_figure_params(dpi=150, figsize=(10, 8))
scv.pl.umap(data, color='cluster_name', legend_loc='none', show=False, alpha=.8)
plt.legend(*plt.gca().get_legend_handles_labels(),
           loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, frameon=False)
plt.subplots_adjust(bottom=0.25)  # adjust spacing to fit legend
plt.savefig(fig_folder + 'X_umap_greenleaf_labels_out.png', bbox_inches='tight')
plt.clf()

scv.settings.set_figure_params(dpi=150, figsize=(8, 6))
scv.pl.umap(data, color='cluster_name', legend_loc='none', show=False)
plt.legend(*plt.gca().get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3)  
plt.tight_layout()
plt.savefig(fig_folder+'X_umap_greenleaf_labels_out.png')
plt.clf()

############################
# Jun.9, plot the sct vector field and reconstructed velocity for greenleaf data
dataset_long = 'greenleaf'
dataset_short = 'glf'

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

split_seed = 317
method = 'sct_GPC'
dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"

adata_prefix = 'adata_'+dataset_short+'_'+method

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

def plot_vf_umap(adata_in,data_version,data,fig_folder,method='sct',celltype_label=None):
    if celltype_label==None: celltype_label = get_celltype_label(data)
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    #scv.tl.velocity_graph(adata)
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    sc.pl.umap(adata, color=celltype_label, ax=axs[0], legend_loc='on data', show=False, frameon=False, legend_fontsize=6)
    print("umapOriginal[0] done")
    sc.pl.umap(adata, color='ptime', ax=axs[1], show=False, frameon=False)
    print("umapOriginal[1] done")
    sct.vf.plot_vector_field(adata,zs_key='X_TNODE',vf_key='X_VF',use_rep_neigh='X_TNODE',color=celltype_label, 
                             show=False,ax=axs[2],legend_loc='none',frameon=False,size=10,alpha=0.4,title=data+' '+data_version,
                             save=fig_folder+'vf/'+data+'_'+method+'_vf_'+data_version+'_umapOriginal_updated.png')
    print("umapOriginal[2] done")    


def plot_metric_withRef(adata,metric,dataset,method,fig_folder,basis_type,split_seed,celltype_label=None,Ngenes=None,recompute=True,basis='umap'):
    metric_color, metric_title = get_metric_color_and_title(metric)
    basis_type = get_basis_type(basis_type)
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    Ngenes_title = ''
    if not Ngenes==None: Ngenes_title = ', Ngenes='+str(Ngenes)
    if basis_type=='Original':
        scv.tl.velocity_graph(adata,n_jobs=8)
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata,basis=basis,color=celltype_label,ax=axs[0],legend_loc='on data',recompute=recompute, legend_fontsize=6,
                                     title="Velocity "+dataset+'+'+method, frameon=False, size=10, alpha=0.5)
    scv.pl.scatter(adata,color=metric_color,cmap='coolwarm',perc=[0,100],ax=axs[1],legend_loc='none', 
                   title=metric_title+" "+dataset+'+'+method+'\n'+Ngenes_title+', split_seed='+str(split_seed),frameon=False,size=10,alpha=0.5)
    plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_withRef_'+basis+basis_type+'_updated.png')
    plt.clf()

def plot_cosine_similarity_withRef(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,recompute=True,celltype_label=None):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    # umapCompute
    adata_plot = adata_total.copy()
    plot_metric_withRef(adata=adata_plot,metric='cos',dataset=dataset,method=method,fig_folder=fig_folder,basis_type='C',
                        split_seed=split_seed,celltype_label=celltype_label,Ngenes=Ngenes,recompute=recompute,basis='umap')
    # umapOriginal
    adata_plot = adata_total.copy()
    adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umapOriginal']
    scv.tl.velocity_graph(adata_plot,n_jobs=8)
    print('Plot umapOriginal')
    plot_metric_withRef(adata=adata_plot,metric='cos',dataset=dataset,method=method,fig_folder=fig_folder,basis_type='Orig',
                        split_seed=split_seed,celltype_label=celltype_label,Ngenes=Ngenes,recompute=recompute,basis='umap')
    

plot_vf_umap(adata_in=total, data_version="total",data=dataset_short,method=method,fig_folder=fig_folder)
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
print('intersected and unioned gene cosine similarity quantiles:\n') # 184, 184

####
# cosine similarity
dataset_long = 'greenleaf'
dataset_short = 'glf'

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

def plot_metric_withRef(adata,metric,dataset,method,fig_folder,basis_type,split_seed,celltype_label=None,Ngenes=None,recompute=True,basis='umap'):
    metric_color, metric_title = get_metric_color_and_title(metric)
    basis_type = get_basis_type(basis_type)
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    Ngenes_title = ''
    if not Ngenes==None: Ngenes_title = ', Ngenes='+str(Ngenes)
    if basis_type=='Original':
        scv.tl.velocity_graph(adata,n_jobs=8)
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata,basis=basis,color=celltype_label,ax=axs[0],legend_loc='on data',recompute=recompute, legend_fontsize=6,
                                     title="Velocity "+dataset+'+'+method, frameon=False, size=10, alpha=0.5)
    scv.pl.scatter(adata,color=metric_color,cmap='coolwarm',perc=[0,100],ax=axs[1],legend_loc='none', 
                   title=metric_title+" "+dataset+'+'+method+'\n'+Ngenes_title+', split_seed='+str(split_seed),frameon=False,size=10,alpha=0.4)
    plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_withRef_'+basis+basis_type+'_updated.png')
    plt.clf()

# sorted
def plot_metric_withRef(adata,metric,dataset,method,fig_folder,basis_type,split_seed,celltype_label=None,Ngenes=None,recompute=True,basis='umap'):
    adata_v0 = adata.copy()
    metric_color, metric_title = get_metric_color_and_title(metric)
    # Get the sort order: largest first, so smallest (blue) last (on top)
    sort_order = adata_v0.obs[metric_color].argsort()#[::-1]
    adata = adata_v0[sort_order].copy() # Reorder everything (obs, X, layers, obsm, etc.)
    # original function steps
    basis_type = get_basis_type(basis_type)
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    Ngenes_title = ''
    if not Ngenes==None: Ngenes_title = ', Ngenes='+str(Ngenes)
    if basis_type=='Original':
        scv.tl.velocity_graph(adata,n_jobs=8)
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata,basis=basis,color=celltype_label,ax=axs[0],legend_loc='on data',recompute=recompute, legend_fontsize=6,
                                     title="Velocity "+dataset+'+'+method, frameon=False, size=10, alpha=0.5)
    scv.pl.scatter(adata,color=metric_color,cmap='coolwarm',perc=[0,100],ax=axs[1],legend_loc='none', 
                   title=metric_title+" "+dataset+'+'+method+'\n'+Ngenes_title+', split_seed='+str(split_seed),frameon=False,size=10,alpha=0.4)
    plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_withRef_'+basis+basis_type+'_updated_sorted.png')
    plt.clf()

# sorted-1
def plot_metric_withRef(adata,metric,dataset,method,fig_folder,basis_type,split_seed,celltype_label=None,Ngenes=None,recompute=True,basis='umap'):
    adata_v0 = adata.copy()
    metric_color, metric_title = get_metric_color_and_title(metric)
    # Get the sort order: largest first, so smallest (blue) last (on top)
    sort_order = adata_v0.obs[metric_color].argsort()#[::-1]
    adata = adata_v0[sort_order].copy() # Reorder everything (obs, X, layers, obsm, etc.)
    # original function steps
    basis_type = get_basis_type(basis_type)
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    Ngenes_title = ''
    if not Ngenes==None: Ngenes_title = ', Ngenes='+str(Ngenes)
    if basis_type=='Original':
        scv.tl.velocity_graph(adata,n_jobs=8)
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata,basis=basis,color=celltype_label,ax=axs[0],legend_loc='on data',recompute=recompute, legend_fontsize=6,
                                     title="Velocity "+dataset+'+'+method, frameon=False, size=10, alpha=0.5)
    scv.pl.scatter(adata,color=metric_color,cmap='coolwarm',perc=[0,100],ax=axs[1],legend_loc='none', 
                   title=metric_title+" "+dataset+'+'+method+'\n'+Ngenes_title+', split_seed='+str(split_seed),frameon=False,size=10,alpha=0.4)
    plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_withRef_'+basis+basis_type+'_updated_sorted-1.png')
    plt.clf()


def plot_cosine_similarity_withRef(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,recompute=True,celltype_label=None):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    # umapCompute
    adata_plot = adata_total.copy()
    plot_metric_withRef(adata=adata_plot,metric='cos',dataset=dataset,method=method,fig_folder=fig_folder,basis_type='C',
                        split_seed=split_seed,celltype_label=celltype_label,Ngenes=Ngenes,recompute=recompute,basis='umap')
    # umapOriginal
    adata_plot = adata_total.copy()
    adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umapOriginal']
    scv.tl.velocity_graph(adata_plot,n_jobs=8)
    print('Plot umapOriginal')
    plot_metric_withRef(adata=adata_plot,metric='cos',dataset=dataset,method=method,fig_folder=fig_folder,basis_type='Orig',
                        split_seed=split_seed,celltype_label=celltype_label,Ngenes=Ngenes,recompute=recompute,basis='umap')

# 1. scv 
split_seed=317
dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'scv_GPC'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

import scvelo as scv
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions_scv import *
from v4_functions import *

celltype_label = 'cluster_name'
total = sc.read(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_total_GPC.h5ad')
split1 = sc.read(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split1_GPC.h5ad')
split2 = sc.read(data_folder+'seed'+str(split_seed)+'/scv_GPC/adata_'+dataset_short+'_'+method+'_split2_GPC.h5ad')
total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
split1.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
split2.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
print('intersected and unioned gene cosine similarity quantiles:\n')

# 2. utv
dataset_long = 'greenleaf'
dataset_short = 'glf'
method = 'utv_GPC'
split_seed=317

import scvelo as scv
import scanpy as sc
#from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
#from v4_functions_scv import *
from v4_functions import *
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_greenleaf/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

total = sc.read(data_folder+'seed'+str(split_seed)+'/utv_GPC/glf_total_GPC_utv.h5ad')
split1 = sc.read(data_folder+'seed'+str(split_seed)+'/utv_GPC/seed317_glf_split1_GPC_utv.h5ad')
split2 = sc.read(data_folder+'seed'+str(split_seed)+'/utv_GPC/seed317_glf_split2_GPC_utv.h5ad')
total.obsm['X_umapOriginal'] = total.obsm['X_umap_greenleaf'].copy()
#compute_umap(total,dataset_short)
#compute_umap(split1,dataset_short)
#compute_umap(split2,dataset_short)

plot_cosine_similarity_withRef(split1,split2,total,dataset_short,method,fig_folder,split_seed)
plot_cosine_similarity_withRef_matpltlib_sorted(adata_split1=split1, adata_split2=split2, adata_total=total, dataset=dataset_short, method=method,
                                                fig_folder=fig_folder, split_seed=split_seed)
plot_cosine_similarity_withRef_matpltlib_sorted_rev(adata_split1=split1, adata_split2=split2, adata_total=total, dataset=dataset_short, method=method,
                                                fig_folder=fig_folder, split_seed=split_seed)


c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
print('intersected and unioned gene cosine similarity quantiles:\n')

# 3. sct
dataset_long = 'greenleaf'
dataset_short = 'glf'

import scanpy as sc
import scvelo as scv
import sctour as sct
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

split_seed = 317
method = 'sct_GPC'
dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+"/"+method+"/"

adata_prefix = 'adata_'+dataset_short+'_'+method

total = sc.read_h5ad(data_folder+adata_prefix+'_total_v4_outputAdded.h5ad') # 
split1 = sc.read_h5ad(data_folder+adata_prefix+'_split1_v4_outputAdded.h5ad') # 
split2 = sc.read_h5ad(data_folder+adata_prefix+'_split2_v4_outputAdded.h5ad') # 

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
print('intersected and unioned gene cosine similarity quantiles:\n')


# 4. velovi
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

split_seed = 317
method = 'velovi_GPC'
dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

split1 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split1_GPC_outputAdded.h5ad')
split2 = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_split2_GPC_outputAdded.h5ad')
total = sc.read_h5ad(data_folder+'adata_glf_velovi_GPC_total_GPC_outputAdded.h5ad')

c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
print('intersected and unioned gene cosine similarity quantiles:\n')

# 5. velovi_woprep
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *

split_seed = 317
method = 'velovi_woprep_GPC'
dataset_short = 'glf'
dataset_long = 'greenleaf'
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_"+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'

split1 = sc.read_h5ad(data_folder+'adata_glf_velovi_woprep_GPC_split1_GPC_outputAdded.h5ad')
split2 = sc.read_h5ad(data_folder+'adata_glf_velovi_woprep_GPC_split2_GPC_outputAdded.h5ad')
total = sc.read_h5ad(data_folder+'adata_glf_velovi_woprep_GPC_total_GPC_outputAdded.h5ad')


c1,n1 = compute_cosine_similarity_intersect(split1,split2,method)
c2,n2 = compute_cosine_similarity_union(split1,split2,method)
print('intersected and unioned gene cosine similarity quantiles:\n')


#####
import matplotlib.pyplot as plt
import numpy as np

def plot_cosine_similarity_withRef_matpltlib_sorted(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,metric_color='cos_sim',celltype_label=None):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    adata = adata_total.copy()
    X = adata.obsm["X_umapOriginal"][:, 0]
    Y = adata.obsm["X_umapOriginal"][:, 1]
    color_values = cos_sim
    sort_idx = np.argsort(color_values)  # Smallest at the end
    X_sorted = X[sort_idx]
    Y_sorted = Y[sort_idx]
    color_sorted = color_values[sort_idx]
    # Plot
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sc = ax.scatter( X_sorted, Y_sorted, c=color_sorted, cmap="coolwarm", vmin=-1, vmax=1, s=4, alpha=0.2, edgecolors="none" )
    ax.set_title(metric_color + " " + dataset + '+' + method + ' (Ngenes=' + str(Ngenes) + ', split_seed=' + str(split_seed)+ ')')
    ax.axis("off")
    plt.colorbar(sc, ax=ax, label='cos_sim')
    plt.tight_layout()
    plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_cos_sim_umapOriginal_sorted_matpltlib.png')
    plt.clf()

def plot_cosine_similarity_withRef_matpltlib_sorted_rev(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,metric_color='cos_sim',celltype_label=None):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    adata = adata_total.copy()
    X = adata.obsm["X_umapOriginal"][:, 0]
    Y = adata.obsm["X_umapOriginal"][:, 1]
    # reverse
    color_values = cos_sim
    sort_idx = np.argsort(color_values)[::-1]  # Smallest at the end
    X_sorted = X[sort_idx]
    Y_sorted = Y[sort_idx]
    color_sorted = color_values[sort_idx]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sc = ax.scatter( X_sorted, Y_sorted, c=color_sorted, cmap="coolwarm", vmin=-1, vmax=1, s=4, alpha=0.2, edgecolors="none" )
    ax.set_title(metric_color + " " + dataset + '+' + method + ' (Ngenes=' + str(Ngenes) + ', split_seed=' + str(split_seed)+ ')')
    ax.axis("off")
    plt.colorbar(sc, ax=ax, label='cos_sim')
    plt.tight_layout()
    plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_cos_sim_umapOriginal_sorted_matpltlib-1.png')
    plt.clf()

plot_cosine_similarity_withRef_matpltlib_sorted(adata_split1=split1, adata_split2=split2, adata_total=total, dataset=dataset_short, method=method,
                                                fig_folder=fig_folder, split_seed=split_seed)
plot_cosine_similarity_withRef_matpltlib_sorted_rev(adata_split1=split1, adata_split2=split2, adata_total=total, dataset=dataset_short, method=method,
                                                fig_folder=fig_folder, split_seed=split_seed)

