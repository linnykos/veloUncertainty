import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import compute_cosine_similarity_union


def compute_df_cos_sim_corr_across_methods(method, dataset_long, dataset_short, type, data_folder, seeds=[317,320,323,326,329]):
    outputAdded = ''
    if ((method=='sct') | ('velovi' in method)): outputAdded = '_outputAdded'
    for i in range(len(seeds)):
        seed = seeds[i]
        s1 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
        s2 = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')
        if ('u' in type): 
            cos_sim = compute_cosine_similarity_union(s1,s2,method)[0]
            df[dataset_short+'_'+method+'_cos_sim_'+str(seed)] = cos_sim
        total = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4'+outputAdded+'.h5ad')
        if (not 'velocity_confidence' in total.obs.columns):
            scv.tl.velocity_confidence(total)
            velo_conf = total.obs['velocity_confidence'].values
            df[dataset_short+'_'+method+'_velo_conf_'+str(seed)] = velo_conf

#total = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(317)+'/'+'scv'+'/adata_'+dataset_short+'_'+'scv'+'_total_v4'+'.h5ad')

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
dataset_long='erythroid'
dataset_short='ery'
df = pd.DataFrame()
for method in ['scv','utv','sct','velovi','velovi_woprep']:
    print(method)
    compute_df_cos_sim_corr_across_methods(method=method, dataset_long=dataset_long, dataset_short=dataset_short, 
                                           type='u', data_folder=data_folder, seeds=[317])

df.to_csv(data_folder+'v4_'+dataset_long+'/'+dataset_short+'_cos_sim_317.csv')

compute_df_cos_sim_corr_across_methods(method='scv', dataset_long=dataset_long, dataset_short=dataset_short, 
                                           type='u', data_folder=data_folder, seeds=[317])
##############################################
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/y2564li/Downloads/proj_scRNA/data/v4_erythroid/ery_cos_sim_317.csv')

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
plt.figure(figsize=(16, 2))
sns.boxplot(x='method', y='value', hue='metric', data=df_melted, palette=custom_palette, flierprops=flierprops,
            order=['utv','sct','velovi.woprep','velovi','scv'])
plt.title('')
plt.xlabel('')
plt.ylabel('Value')
plt.legend()
plt.xticks()
plt.tight_layout()
ax = plt.gca()  # Get the current axis
for spine in ax.spines.values():
    spine.set_edgecolor('lightgray') 
    spine.set_linewidth(1.5) 

plt.tight_layout()
plt.savefig('/Users/y2564li/Downloads/proj_scRNA/data/v4_erythroid/boxplot.png')

plt.show()