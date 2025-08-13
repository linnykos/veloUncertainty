import scanpy as sc
import pandas as pd
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from v4_functions_transMat import *
import os
import re

dataset_short = 'glf'
dataset_long = 'greenleaf'
celltype_label = 'cluster_name'

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'

## GPC
gene_set_name_prefix = 'GPC'
gene_set_name = gene_set_name_prefix
split_seed = 317


def get_total_split12(method_prefix, gene_set_name):
    method = method_prefix + '_' + gene_set_name
    print(method)
    adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
    adata_path1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('split1', s)][0]
    adata_path2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('split2', s)][0]
    adata_total_path = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('total', s)][0]
    adata_split1 = sc.read_h5ad( adata_path1 )
    adata_split2 = sc.read_h5ad( adata_path2 )
    adata_total = sc.read_h5ad( adata_total_path )
    return adata_total, adata_split1, adata_split2

scv_total, scv_split1, scv_split2 = get_total_split12(method_prefix='scv', gene_set_name=gene_set_name_prefix)
velovi_woprep_total, velovi_woprep_split1, velovi_woprep_split2 = get_total_split12(method_prefix='velovi_woprep', gene_set_name=gene_set_name_prefix)

"""
from numpy.linalg import norm
v_scv = scv_total.layers['velocity']
v_scv = v_scv[:,np.where(~np.isnan(v_scv[0]))[0] ]
v_velovi_woprep = velovi_woprep_total.layers['velocity']
"""

def compute_intersected_genes(adata_split1, adata_split2, method1, method2):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    if "scv" in method1:
        print('scv in method1')
        velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
    if "scv" in method2:
        print('scv in method2')
        velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Number of overlapped genes for velocity computation in two methods = '+str(common_genes_velocity.shape[0])) 
    velo_df1 = pd.DataFrame(adata_split1.layers['velocity'], columns=adata_split1.var.index.tolist())
    velo_df2 = pd.DataFrame(adata_split2.layers['velocity'], columns=adata_split2.var.index.tolist())
    return velo_df1[common_genes_velocity], velo_df2[common_genes_velocity] 

velo_df1, velo_df2 = compute_intersected_genes(adata_split1=scv_total, adata_split2=velovi_woprep_total, method1='scv', method2='v')

# compute angles
from numpy.linalg import norm
angles = np.degrees(np.arccos((velo_df1*velo_df2).sum(1) / (norm(velo_df1, axis=1)*norm(velo_df2, axis=1)) ))

"""
np.mean(angles < 90)
np.where(angles < 45)

np.mean(angles > 90)
np.where(angles > 45)

scv_total[np.where(angles < 45)[0],:].obs['cluster_name']
scv_total[np.where(angles > 90)[0],:].obs['cluster_name']
"""

# proportion in each celltype, angle>90
from collections import Counter
for i in scv_total.obs['cluster_name'].cat.categories:
    prop_gt90 = np.round(Counter(scv_total[np.where(angles > 90)[0],:].obs['cluster_name'])[i]/ Counter(scv_total.obs['cluster_name'])[i], 4)
    print(i+": "+ str(prop_gt90))

"""
Counter(scv_total[np.where(angles > 90)[0],:].obs['cluster_name'])
Counter(scv_total.obs['cluster_name'])

[i for i in scv_total.obs['cluster_name'].cat.categories]
[Counter(scv_total[np.where(angles < 90)[0],:].obs['cluster_name'])[i]/ Counter(scv_total.obs['cluster_name'])[i] for i in scv_total.obs['cluster_name'].cat.categories]
"""

scv_total.obs['angle_large'] = np.array(angles>90)
scv_total.obs['velo_angle'] = np.array(angles)

adata1 = scv_total.copy()
del adata1.obsm['X_umap']

# plot (angle>90)
scv.pl.scatter(adata1, color='angle_large', basis='umap_greenleaf',title='velo angle scv vs velovi_woprep, glf GPC', perc=[0, 100], 
               save=fig_folder+'velo_angle_scv_vs_velovi_woprep_scatter.png', size=5)
scv.pl.scatter(adata1, color='angle_large', basis='X_umap',title='velo angle scv vs velovi_woprep, glf GPC', perc=[0, 100], 
               save=fig_folder+'velo_angle_scv_vs_velovi_woprep_scatter_umapCompute.png', size=5)
# plot angle, umap original
scv.pl.scatter(adata1[adata1.obs.sort_values('velo_angle').index], color='velo_angle', basis='umap_greenleaf',
               title='velo angle scv vs velovi_woprep, glf GPC', perc=[0, 100], cmap='coolwarm',
               save=fig_folder+'velo_angle_num_scv_vs_velovi_woprep_scatter.png', size=5)

# plot umap of clusters
scv.pl.scatter(adata1, color='cluster_name', basis='umap_greenleaf',title='scv umap, glf GPC', legend_loc='right margin',
               save=fig_folder+'scv_scatter.png', size=5)


import matplotlib.pyplot as plt
cell_cat = scv_total.obs[celltype_label].cat.categories
fig, ax = plt.subplots(figsize=(len(cell_cat)*1.6, 9))
#data_to_plot = [scv_total.obs['velo_angle'][scv_total.obs[celltype_label].array==celltype] for celltype in cell_cat]

data_to_plot = [angles[scv_total.obs[celltype_label].array==celltype] for celltype in cell_cat]

# Create boxplot
ax.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'),
            medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
            flierprops=dict(marker='o', color='rosybrown', alpha=0.5, markerfacecolor='aliceblue', markeredgecolor='rosybrown'))
counts = [len(data) for data in data_to_plot]
x_labels = [f'{cell_cat[i]} (n={counts[i]})' for i in range(len(cell_cat))]
ax.set_xticks(range(1, len(cell_cat) + 1))
ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=12)
ax.set_xlabel('Cell Types')
ax.set_ylabel('velocity angle')
ax.set_ylim(0,180)
ax.set_title(f'velo angle scv vs velovi_woprep, glf GPC')
plt.tight_layout()
plt.savefig(fig_folder+'velo_angle_num_scv_vs_velovi_woprep_byCelltype_boxplot1.png')
plt.close()


## histogram
plt.clf()
cell_cat = scv_total.obs[celltype_label].cat.categories
fig,axs = plt.subplots(ncols=4, nrows=3, figsize=(25,16))
axs = axs.ravel()
for idx,ax in enumerate(axs):
    if idx==len(cell_cat): break
    celltype = cell_cat[idx]
    cos_sim_celltype = angles[scv_total.obs[celltype_label].array==celltype]
    Ncells = cos_sim_celltype.shape[0]
    counts, bins, patches = ax.hist(cos_sim_celltype, bins=20, edgecolor='gainsboro',color='powderblue') 
    max_frequency = np.max(counts)
    ax.axvline(np.mean(cos_sim_celltype), color='brown', linestyle='dashed', linewidth=1.5)
    ax.axvline(np.median(cos_sim_celltype), color='peru', linestyle='dashed', linewidth=1.5)
    text_x = np.quantile(cos_sim_celltype,[.0])[0]
    text_y = max_frequency/5
    ax.text(text_x,text_y*3,'mean='+str(np.round(np.mean(cos_sim_celltype),4)), color='firebrick', fontsize=11)
    ax.text(text_x,text_y*2,'median='+str(np.round(np.median(cos_sim_celltype),4)), color='sienna', fontsize=11)
    ax.set_xlabel('velo angle')
    ax.set_ylabel('Frequency')
    ax.set_title(celltype+' (Ncells='+str(Ncells))

plt.savefig(fig_folder+'velo_angle_num_scv_vs_velovi_woprep_byCelltype_hist.png')
plt.clf()



########
## nGPCblk227
gene_set_name_prefix = 'nGPCblk'
grid_seed = 227
gene_set_name = gene_set_name_prefix + str(grid_seed)
split_seed = 317

# nGPCblk230
gene_set_name_prefix = 'nGPCblk'
grid_seed = 230
gene_set_name = gene_set_name_prefix + str(grid_seed)
split_seed = 320

def get_total_split12(method_prefix, gene_set_name):
    method = method_prefix + '_' + gene_set_name
    print(method)
    adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
    if 'scv' in method:
        adata_path1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('split1.*h5ad', s)][0]
        adata_path2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('split2.*h5ad', s)][0]
        adata_total_path = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('total.*h5ad', s)][0]
    else:
        adata_path1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('split1.*outputAdded.*h5ad', s)][0]
        adata_path2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('split2.*outputAdded.*h5ad', s)][0]
        adata_total_path = data_folder+'seed'+str(split_seed)+'/'+method+'/'+[s for s in adata_names if re.search('total.*outputAdded.*h5ad', s)][0]
    adata_split1 = sc.read_h5ad( adata_path1 )
    adata_split2 = sc.read_h5ad( adata_path2 )
    adata_total = sc.read_h5ad( adata_total_path )
    return adata_total, adata_split1, adata_split2

scv_total, scv_split1, scv_split2 = get_total_split12(method_prefix='scv', gene_set_name=gene_set_name)
velovi_woprep_total, velovi_woprep_split1, velovi_woprep_split2 = get_total_split12(method_prefix='velovi_woprep', gene_set_name=gene_set_name)

velo_df1, velo_df2 = compute_intersected_genes(adata_split1=scv_total, adata_split2=velovi_woprep_total, method1='scv', method2='v')

from numpy.linalg import norm
angles = np.degrees(np.arccos((velo_df1*velo_df2).sum(1) / (norm(velo_df1, axis=1)*norm(velo_df2, axis=1)) ))

from collections import Counter
for i in scv_total.obs['cluster_name'].cat.categories:
    prop_gt90 = np.round(Counter(scv_total[np.where(angles > 90)[0],:].obs['cluster_name'])[i]/ Counter(scv_total.obs['cluster_name'])[i], 4)
    print(i+": "+ str(prop_gt90))

scv_total.obs['angle_large'] = np.array(angles>90)
scv_total.obs['velo_angle'] = np.array(angles)

adata1 = scv_total.copy()
del adata1.obsm['X_umap']

# plot (angle>90)
scv.pl.scatter(adata1, color='angle_large', basis='umap_greenleaf',title='velo angle scv vs velovi_woprep, glf nGPCblk227', perc=[0, 100], 
               save=fig_folder+'velo_angle_nGPC_scv_vs_velovi_woprep_scatter.png', size=5)
# plot angle, umap original
scv.pl.scatter(adata1[adata1.obs.sort_values('velo_angle').index], color='velo_angle', basis='umap_greenleaf',
               title='velo angle scv vs velovi_woprep, glf nGPCblk227', perc=[0, 100], cmap='coolwarm',
               save=fig_folder+'velo_angle_nGPC_num_scv_vs_velovi_woprep_scatter.png', size=5)

## nGPCblk230
# plot (angle>90)
scv.pl.scatter(adata1[adata1.obs.sort_values('angle_large').index], color='angle_large', basis='umap_greenleaf',title='velo angle scv vs velovi_woprep, glf nGPCblk230', perc=[0, 100], 
               save=fig_folder+'velo_angle_nGPCblk230_scv_vs_velovi_woprep_scatter.png', size=5)
# plot angle, umap original
scv.pl.scatter(adata1[adata1.obs.sort_values('velo_angle').index], color='velo_angle', basis='umap_greenleaf',
               title='velo angle scv vs velovi_woprep, glf nGPCblk230', perc=[0, 100], cmap='coolwarm',
               save=fig_folder+'velo_angle_nGPCblk230_num_scv_vs_velovi_woprep_scatter.png', size=5)

