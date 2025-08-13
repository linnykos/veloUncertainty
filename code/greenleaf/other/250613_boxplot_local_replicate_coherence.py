import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import compute_cosine_similarity_union
import os
import re

def compute_df_cos_sim_corr_across_methods(method, dataset_long, dataset_short, type, data_folder, seeds=[317,320,323,326,329]):
    outputAdded = ''
    if ((method=='sct') | ('velovi' in method)): outputAdded = '_outputAdded'
    for i in range(len(seeds)):
        split_seed = seeds[i]
        print(split_seed)
        adata_names = os.listdir(data_folder+'seed'+str(split_seed)+'/'+method)
        if 'scv' in method or 'utv' in method:
            path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1', s)][0]
            path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2', s)][0]
            path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total', s)][0]
        else: 
            path_split1 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split1.*outputAdded', s)][0]
            path_split2 = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('split2.*outputAdded', s)][0]
            path_total = data_folder+'seed'+str(split_seed)+'/'+method+'/'+ [s for s in adata_names if re.search('total.*outputAdded', s)][0]
        s1 = sc.read_h5ad( path_split1 )
        s2 = sc.read_h5ad( path_split2 )
        if ('u' in type): 
            cos_sim = compute_cosine_similarity_union(s1,s2,method)[0]
            df[method+'_replicate_coherence_'+str(split_seed)] = cos_sim
        total = sc.read_h5ad( path_total )
        if (not 'velocity_confidence' in total.obs.columns):
            scv.tl.velocity_confidence(total)
            velo_conf = total.obs['velocity_confidence'].values
            df[method+'_local_coherence_'+str(split_seed)] = velo_conf

dataset_long = 'greenleaf'
dataset_short = 'glf'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

df = pd.DataFrame()
for method in ['scv','utv','sct','velovi','velovi_woprep']:
    print(method)
    compute_df_cos_sim_corr_across_methods(method=method, dataset_long=dataset_long, dataset_short=dataset_short, 
                                           type='u', data_folder=data_folder, seeds=[317])

df.to_csv(data_folder+dataset_short+'_local_rep_coherence_317.csv')

adata = sc.read_h5ad(data_folder+'/glf_total_allgenes.h5ad')
spliced = adata.layers['spliced'].A
spliced[spliced != 0].mean() # Average non-zero spliced counts across cells
np.mean( np.array(spliced.sum(axis=0)).flatten()>0 ) # Fraction of genes with non-zero spliced expression


############################
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

dataset_long='greenleaf'
dataset_short='glf'
fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/'
data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'

df = pd.read_csv(data_folder+dataset_short+'_local_rep_coherence_317.csv')

df.columns = ['Unnamed: 0', 'scv_replicate_coherence', 'scv_local_coherence',
       'utv_replicate_coherence', 'utv_local_coherence', 'sct_replicate_coherence',
       'sct_local_coherence', 'velovi_replicate_coherence',
       'velovi_local_coherence', 'velovi.woprep_replicate_coherence',
       'velovi.woprep_local_coherence']

custom_palette = {'replicate_coherence': '#ebb5ae', 'local_coherence': '#fff8a8'}  
df_melted = pd.melt(df, id_vars='Unnamed: 0', var_name='method_metric', value_name='value')
df_melted['method'] = df_melted['method_metric'].apply(lambda x: '_'.join(x.split('_')[:1]))
df_melted['metric'] = df_melted['method_metric'].apply(lambda x: '_'.join(x.split('_')[1:]))
flierprops = dict(marker='o', markersize=3, linestyle='none', markerfacecolor='none', markeredgecolor=(0.7, 0.7, 0.7, 0.05)) 
plt.figure(figsize=(16, 4))
sns.boxplot(x='method', y='value', hue='metric', data=df_melted, palette=custom_palette, flierprops=flierprops,
            order=['utv','sct','velovi.woprep','velovi','scv'])

plt.legend(loc='lower right')
plt.title('Boxplots of local vs replicate coherence ('+dataset_long+')', fontsize=16)
plt.xlabel('')
plt.ylabel('Value', fontsize=14)
plt.legend()
plt.xticks(fontsize=14)
plt.tight_layout()
ax = plt.gca()  # Get the current axis
for spine in ax.spines.values():
    spine.set_edgecolor('lightgray') 
    spine.set_linewidth(1.5) 

plt.tight_layout()
plt.savefig(fig_folder+'local_rep_coherence_boxplot.png')
plt.clf()

