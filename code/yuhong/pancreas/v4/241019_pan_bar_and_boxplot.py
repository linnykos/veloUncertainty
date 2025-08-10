import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# pancreas
data = {
    'Method': ['utv', 'sct', 'velovi.woprep', 'velovi', 'scv'],
    'Positive & significant': ['10/10', '10/10', '6/10', '10/10', '8/10']
}

df = pd.DataFrame(data)
df['Positive'] = df['Positive & significant'].str.split('/').str[0].astype(int)
plt.figure(figsize=(6, 4.5))
sns.barplot(x='Method', y='Positive', data=df, color='#a8c9ff')
plt.title('')
plt.xlabel('')
plt.ylabel('Frequency')
plt.xticks()
ax = plt.gca()  # Get the current axis
for spine in ax.spines.values():
    spine.set_edgecolor('lightgray') 
    spine.set_linewidth(1.5) 

plt.tight_layout()
plt.savefig('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/bar_cos_sim_corr_pval.png')
plt.show()

# pancreasINC
data = {
    'Method': ['utv', 'sct', 'velovi.woprep', 'velovi', 'scv'],
    'Positive & significant': ['6/10', '6/10', '7/10', '10/10', '8/10']
}
df = pd.DataFrame(data)
df['Positive'] = df['Positive & significant'].str.split('/').str[0].astype(int)
plt.figure(figsize=(6, 4.5))
sns.barplot(x='Method', y='Positive', data=df, color='#a8c9ff')
plt.title('')
plt.xlabel('')
plt.ylabel('Frequency')
plt.xticks()
ax = plt.gca()  # Get the current axis
for spine in ax.spines.values():
    spine.set_edgecolor('lightgray') 
    spine.set_linewidth(1.5) 

plt.tight_layout()
plt.savefig('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreasINC/bar_cos_sim_corr_pval.png')


###############################
import scanpy as sc
## pan
df = pd.read_csv('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/pan_cos_sim_corr.csv')
data_to_plot = df[['pan_utv_317','pan_sct_317','pan_velovi_woprep_317','pan_velovi_317','pan_scv_317']]

plt.figure(figsize=(6, 5))
plt.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'), 
            medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
            flierprops=dict(marker='o', markersize=3, linestyle='none', color='rosybrown', alpha=0.05, 
                            markerfacecolor='none', markeredgecolor=(0.7, 0.7, 0.7, 0.05)))
plt.title('')
plt.xticks([1, 2, 3, 4, 5], ['utv', 'sct', 'velovi_woprep','velovi', 'scv'])
plt.ylim(-1, 1)
plt.savefig('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/boxplot.png')

## panINC
df = pd.read_csv('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreasINC/panINC_cos_sim_corr.csv')
data_to_plot = df[['panINC_utv_320','panINC_sct_320','panINC_velovi_woprep_320','panINC_velovi_320','panINC_scv_320']]

plt.figure(figsize=(6, 5))
plt.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'), 
            medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
            flierprops=dict(marker='o', markersize=3, linestyle='none', color='rosybrown', alpha=0.5, 
                            markerfacecolor='none', markeredgecolor=(0.7, 0.7, 0.7, 0.05)))
plt.title('')
plt.xticks([1, 2, 3, 4, 5], ['utv', 'sct', 'velovi_woprep','velovi', 'scv'])
plt.ylim(-1, 1)
plt.savefig('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreasINC/boxplot.png')

################
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_scatter_adjustColor(adata,metric,plot_type,fig_folder,dataset,method,basis_type,split_seed,basis='umap',Ngenes=None,recompute=True,cmap='coolwarm'):
    metric_color, metric_title = get_metric_color_and_title(metric)
    basis_type = get_basis_type(basis_type)
    Ngenes_title = ''
    if not Ngenes==None: Ngenes_title = ', Ngenes='+str(Ngenes)
    if basis_type=='Original': scv.tl.velocity_graph(adata,n_jobs=8)
    if 'emb' in plot_type:
        scv.pl.velocity_embedding_stream(adata, basis=basis, color=metric_color, cmap=cmap, recompute=recompute,
                                         perc=[1, 100], title=metric_title+" "+dataset+'+'+method+'\n'+Ngenes_title+', split_seed='+str(split_seed),
                                         save=fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_'+basis+basis_type+'.png')
    elif 'sca' in plot_type:
        norm = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0.25, vmax=1)
        cmap = mpl.cm.coolwarm
        scv.pl.scatter(adata, color=metric_color, cmap=cmap, title=metric_title + " " + dataset + '+' + method + '\n' + Ngenes_title + ', split_seed=' + str(split_seed), 
                    save=fig_folder + "metric/" + dataset + '_' + method + '_' + metric_color + '_' + basis + basis_type + '_scatter_adjustColor.png',
                    perc=[0, 100], norm=norm)

plot_scatter_adjustColor()

data_folder = '/Users/y2564li/Downloads/proj_scRNA/data/'
dataset_long = 'pancreas'
dataset_short = 'pan'
seed = 317
method = 'scv'
total = sc.read_h5ad(data_folder+'v4_'+dataset_long+'/seed'+str(seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_v4'+'.h5ad')

df = pd.read_csv('/Users/y2564li/Downloads/proj_scRNA/data/v4_pancreas/pan_cos_sim_corr.csv')
total.obs['cos_sim'] = df['pan_scv_317'].values
norm = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0.25, vmax=1)
total.obs['cos_sim_norm'] = norm(total.obs['cos_sim'])
cmap = mpl.cm.coolwarm
sc.pl.scatter(total, color='cos_sim_norm', color_map=cmap, basis='umap')


