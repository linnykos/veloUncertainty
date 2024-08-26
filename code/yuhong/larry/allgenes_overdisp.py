import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import scvelo as scv
import statsmodels.formula.api as smf
import pandas as pd

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v2_functions import compute_gene_correlation_between_splits
from countsplit import estimate_overdisps

"""
# test - random vector
np.random.seed(42)
nb_sample = np.random.negative_binomial(n=1, p=.01, size=1000)
estimate_overdisps(np.matrix(nb_sample)) # 0.98728751
"""

# Estimate the overdispersion using the negative binomial model
"""
# test - ChatGPT ver overdispersion
result = None
theta = None
overdispersion = None
nb_sample = S[:,3].todense().flatten()
data = pd.DataFrame({'y': nb_sample.A1})
model = smf.negativebinomial('y ~ 1', data)
result = model.fit()
theta = result.params[-1]
overdispersion = 1 / theta
overdispersion
"""

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
df = pd.read_csv(data_folder+'v2_larry/larry_df_allgenes.csv')
#total2 = ad.read_h5ad(data_folder+"v2_larry/larry.h5ad")
#S = total2.layers['spliced']
#estimate_overdisps(S[:,0:15])
# array([0.95122942,        inf,        inf, 0.53641167,        inf])

"""
# test - manually implemented overdispersion estimation function
y = S[:,0]
if hasattr(y, 'todense'): y = y.todense().A 
df = pd.DataFrame({'counts': y.flatten()})
model = smf.negativebinomial('counts ~ 1', data=df)
result = model.fit()
alpha = result.params['alpha']
b = 1/alpha
b

# test - random Poisson sample
np.random.seed(2410)
samples = np.random.poisson(lam=5, size=10000).reshape(2000,5)
estimate_overdisps(samples)
# array([6.54010220e+05, 2.76667314e+05, 1.29765115e+06, 1.25197669e+02, 3.80863686e+05])
"""

# sct genes
dataset_long = 'larry'
dataset_short = 'larry'
method = 'sct'
total=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_total_v2.h5ad')
split1=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_split1_v2.h5ad')
split2=sc.read_h5ad(data_folder+'v2_'+dataset_long+'/'+method+'/adata_newvelo_'+dataset_short+'_'+method+'_split2_v2.h5ad')
adata = sc.read_h5ad(data_folder+"v2_larry/larry.h5ad") # n_obs × n_vars = 49302 × 23420

import matplotlib as mpl
def rgb2hex(rgb):
    r = int(rgb[0]*255)
    g = int(rgb[1]*255)
    b = int(rgb[2]*255)
    return "#{:02x}{:02x}{:02x}".format(r,g,b)

split1.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[0:415:40]]
split2.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[0:415:40]]
total.uns['state_info_colors'] = [rgb2hex(color) for color in mpl.colormaps['twilight_shifted'].colors[0:415:40]]

"""
gene_names_adata = adata.var.index.copy()
positions_dict_adata = {gene: pos for pos, gene in enumerate(gene_names_adata)}
df['gene_names']=gene_names_adata

def compute_nnz(adata,type):
    n_rows = adata.layers[type].shape[0]
    nonzeros_per_column = adata.layers[type].getnnz(axis=0)
    fraction_nonzeros = nonzeros_per_column / n_rows
    return fraction_nonzeros

df['frac_nnz_s'] = compute_nnz(adata,'spliced')
df['frac_nnz_u'] = compute_nnz(adata,'unspliced')
"""
common_genes = np.intersect1d(np.array(split1.var.index), np.array(split2.var.index))
df_sct = df.loc[df['gene_names'].isin(common_genes)]

# filter out genes that are highly-correlated and have large overdispersion estimates
np.intersect1d(df_sct['gene_names'][df_sct['overdisps_S']>10],df_sct['gene_names'][df_sct['cor_spliced']>.8])
#array(['Akr1c18', 'BC100530', 'H2-Aa', 'Itgb3', 'Nrgn', 'Parvb', 'Ppbp', 'Saa3', 'Thbs1', 'Timp3'], dtype=object)

fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v2_larry/" 
import matplotlib.patches as patches
def plot_gene_exp2splits(split1,split2,gene_name):
    gene_idx1 = np.where(split1.var.index==gene_name)[0][0]
    gene_idx2 = np.where(split2.var.index==gene_name)[0][0]
    color_map = dict(zip(total.obs['state_info'].cat.categories, total.uns['state_info_colors'])) 
    colors = [color_map[i] for i in total.obs['state_info']]
    x_s = split1.layers['spliced'][:,gene_idx1].todense().A1
    y_s = split2.layers['spliced'][:,gene_idx2].todense().A1
    x_u = split1.layers['unspliced'][:,gene_idx1].todense().A1
    y_u = split2.layers['unspliced'][:,gene_idx2].todense().A1
    np.random.seed(1533)
    x_s = x_s + np.clip(np.random.normal(loc=0, scale=1, size=len(x_s)), -1, 1)
    y_s = y_s + np.clip(np.random.normal(loc=0, scale=1, size=len(y_s)), -1, 1)
    x_u = x_u + np.clip(np.random.normal(loc=0, scale=1, size=len(x_u)), -1, 1)
    y_u = y_u + np.clip(np.random.normal(loc=0, scale=1, size=len(y_u)), -1, 1)
    corr_s = np.round(df['cor_spliced'][df['gene_names']==gene_name].values[0],3)
    corr_u = np.round(df['cor_unspliced'][df['gene_names']==gene_name].values[0],3)
    overdisp_s = df['overdisps_S'][df['gene_names']==gene_name].values[0]
    overdisp_u = df['overdisps_U'][df['gene_names']==gene_name].values[0]
    nnz_s = np.round(df['frac_nnz_s'][df['gene_names']==gene_name].values[0],3)
    nnz_u = np.round(df['frac_nnz_u'][df['gene_names']==gene_name].values[0],3)
    handles = [patches.Patch(color=color, label=label) for label, color in color_map.items()]
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
    axes[0].scatter(x_s, y_s, color=colors, alpha=0.2)
    axes[0].legend(handles=handles, bbox_to_anchor=(1.1, 0), loc='lower right', framealpha=0.3)
    axes[0].patch.set_alpha(0.5)
    axes[0].set_title('Scatter plot of gene expression in splits ('+gene_name+', spliced)\n'+
                        'corr='+str(corr_s)+', overdisp='+f"{overdisp_s:.3e}"+', nonzero%='+str(nnz_s))
    axes[0].set_xlabel("split1")
    axes[0].set_ylabel("split2")
    axes[1].scatter(x_u, y_u, color=colors, alpha=0.2)
    axes[1].legend(handles=handles, bbox_to_anchor=(1.1, 0), loc='lower right', framealpha=0.3)
    axes[1].patch.set_alpha(0.5)
    axes[1].set_title('Scatter plot of gene expression in splits ('+gene_name+', unspliced)\n'+
                        'corr='+str(corr_u)+', overdisp='+f"{overdisp_u:.3e}"+', nonzero%='+str(nnz_u))
    axes[1].set_xlabel("split1")
    axes[1].set_ylabel("split2")
    plt.tight_layout()
    plt.savefig(fig_folder+'gene_'+gene_name+'.png') 
    plt.clf()
    print(gene_name+' done!')

"""
plt.scatter(x_s, y_s, color=colors, alpha=0.2)
plt.legend(handles=handles, bbox_to_anchor=(1.1, 0), loc='lower right', framealpha=0.3)
plt.title('Scatter plot of gene expression in splits ('+gene_name+', spliced)\n'+'corr='+str(corr)+', overdisp='+f"{overdisp:.3e}"+', nonzero%='+str(nnz))
plt.xlabel('split1')
plt.ylabel('split2')
plt.savefig(fig_folder+'gene_'+gene_name+'.png') 
plt.clf()
"""

plot_gene_exp2splits(split1,split2,'Akr1c18')

for gene_name in np.intersect1d(df_sct['gene_names'][df_sct['overdisps_S']>10],df_sct['gene_names'][df_sct['cor_spliced']>.8]):
    plot_gene_exp2splits(split1,split2,gene_name)

gene_names_intersect = np.intersect1d(df_sct['gene_names'][df_sct['overdisps_S']>10],
                                      df_sct['gene_names'][df_sct['cor_spliced']>.8])

np.intersect1d(df['gene_names'][df['overdisps_S']>10],
               df['gene_names'][df['cor_spliced']>.8])

[np.where(df['gene_names']==gene_name)[0][0] for gene_name in gene_names_intersect]
# [4946, 7739, 8544, 4005, 21172, 7411, 16169, 18303, 12291, 2191]
overdisp_smf = []
overdisp_countsplit = []
overdisp_df = []
for gene_name in gene_names_intersect:
    gene_idx = np.where(df['gene_names']==gene_name)[0][0]
    print(gene_name, gene_idx, np.sum(adata.layers['spliced'][:,gene_idx].todense().A1))
    result = None
    theta = None
    overdispersion = None
    nb_sample = adata.layers['spliced'][:,gene_idx].todense().flatten()
    data = pd.DataFrame({'y': nb_sample.A1})
    model = smf.negativebinomial('y ~ 1', data)
    result = model.fit()
    theta = result.params[-1]
    overdispersion = 1 / theta
    overdisp_smf.append(np.round(overdispersion,5))
    overdisp_df.append(np.round(df['overdisps_S'][gene_idx],5))
    overdisp_countsplit.append(np.round(estimate_overdisps(adata.layers['spliced'][:,gene_idx]),5))

df['overdisps_S'][np.where(df['gene_names'] in gene_names_intersect)]

idx = 7739
idx_start = idx-2
idx_end = idx+3
estimate_overdisps(adata.layers['spliced'][:,idx_start:idx_end])
df['overdisps_S'][idx_start:idx_end]

estimate_overdisps(adata.layers['spliced'][:,idx])

estimate_overdisps(adata.layers['spliced'][:,7739:7750])
estimate_overdisps(adata.layers['spliced'][:,7737:7740])
estimate_overdisps(adata.layers['spliced'][:,7720:7740])
df['overdisps_S'][7739]
"""
# test
total2.var.index
total2.layers['spliced']
gene_name = 'Saa3'
idx = np.where(total2.var.index == gene_name)[0][0]
arr = total2.layers['spliced'][:,idx].todense().flatten()
print('fraction of nonzeros = ',np.sum(arr>0)/total2.shape[0])
print('mean across nonzero entries', np.mean(arr[arr>0]))
print('var across nonzero entries', np.var(arr[arr>0]))

arr[np.where(arr>0)].tolist() # nonzero elements
"""

np.corrcoef(split1.layers['spliced'][:,np.where(split1.var.index=='Timp3')[0][0]].todense().A1,
            split2.layers['spliced'][:,np.where(split2.var.index=='Timp3')[0][0]].todense().A1)
df['cor_spliced'][2191]

np.corrcoef(split1.layers['spliced'][:,np.where(split1.var.index=='Saa3')[0][0]].todense().A1,
            split2.layers['spliced'][:,np.where(split2.var.index=='Saa3')[0][0]].todense().A1)
df['cor_spliced'][18303]

def compute_corr_and_compare(split1,split2,gene_name,splice='spliced'):
    corr = np.corrcoef(split1.layers[splice][:,np.where(split1.var.index==gene_name)[0][0]].todense().A1,
                        split2.layers[splice][:,np.where(split2.var.index==gene_name)[0][0]].todense().A1)
    corr_df = df['cor_'+splice][df['gene_names']==gene_name].values[0]
    return np.round(corr[0,1],6)-np.round(corr_df,6)

compute_corr_and_compare(split1,split2,common_genes[0])
[compute_corr_and_compare(split1,split2,common_genes[i]) for i in range(50)]

from countsplit import countsplit
np.random.seed(317)
tmp1,tmp2 = countsplit(X=adata.layers['spliced'][:,7739], folds=2, epsilon=None, overdisps=np.array([0.01472]))

color_map = dict(zip(total.obs['state_info'].cat.categories, total.uns['state_info_colors'])) 
colors = [color_map[i] for i in total.obs['state_info']]

import matplotlib.patches as patches
def plot_gene_exp2splits_nb(gene_name,overdisp):
    gene_idx = np.where(adata.var.index==gene_name)[0][0]
    np.random.seed(317)
    tmp1,tmp2 = countsplit(X=adata.layers['spliced'][:,gene_idx], folds=2, epsilon=None, overdisps=np.array([overdisp]))
    plt.scatter(tmp1.todense().A1, tmp2.todense().A1, color=colors, alpha=0.2)
    handles = [patches.Patch(color=color, label=label) for label, color in color_map.items()]
    plt.legend(handles=handles, bbox_to_anchor=(1.1, 0), loc='lower right', framealpha=0.3)
    #plt.patch.set_alpha(0.5)
    plt.title('Scatter plot of gene expression in splits ('+gene_name+', spliced)\n'+
              'overdisp='+f"{overdisp:.3e}"+', correlation='+str(np.round(np.corrcoef(tmp1.todense().A1,tmp2.todense().A1)[0,1],4)))
    plt.xlabel("split1 (new overdisp)")
    plt.ylabel("split2 (new overdisp)")
    plt.savefig(fig_folder+'gene_'+gene_name+'_nb.png')
    plt.clf()


plot_gene_exp2splits_nb(gene_name='Akr1c18',overdisp=0.02677)
plot_gene_exp2splits_nb(gene_name='BC100530',overdisp=0.01472)
plot_gene_exp2splits_nb(gene_name='H2-Aa',overdisp=0.02338)
plot_gene_exp2splits_nb(gene_name='Itgb3',overdisp=0.11167)
plot_gene_exp2splits_nb(gene_name='Nrgn',overdisp=0.04765)
plot_gene_exp2splits_nb(gene_name='Parvb',overdisp=0.02950)
plot_gene_exp2splits_nb(gene_name='Ppbp',overdisp=0.01481)
plot_gene_exp2splits_nb(gene_name='Saa3',overdisp=0.06224)
plot_gene_exp2splits_nb(gene_name='Thbs1',overdisp=0.02089)
plot_gene_exp2splits_nb(gene_name='Timp3',overdisp=0.01280)



