import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
import datetime

## print current time with a message
def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

### used in 2splitCorr.py, helper function inside plot_gene_correlation_between_splits
def compute_gene_correlation_between_splits(mat1,mat2):
    if hasattr(mat1, 'todense'):
        mat1 = mat1.todense().A 
        mat2 = mat2.todense().A 
    m1 = np.transpose(mat1)
    m2 = np.transpose(mat2)
    cor = []
    for i in range(m1.shape[0]):
        if i%1000==0:
            print(i)
        if np.sum(m1[i,:])==0 and np.sum(m2[i,:])==0:
            cor.append(np.nan)
        else:
            cor.append(np.corrcoef(m1[i,:],m2[i,:])[0,1])
    cor = np.array(cor)
    print("Number of valid values = "+str(cor[~np.isnan(cor)].shape[0]))
    print("Quantiles: "+str(np.quantile(cor[~np.isnan(cor)],[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])))
    return cor

### used in 2splitCorr.py, to plot correlations between splits
def plot_gene_correlation_between_splits(adata1,adata2,fig_path,fig_folder):
    common_genes = np.intersect1d(np.array(adata1.var.index), np.array(adata2.var.index))
    gene_names_split1 = adata1.var.index.copy()
    positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
    positions_split1 = [positions_dict_split1[gene] for gene in common_genes]
    gene_names_split2 = adata2.var.index.copy()
    positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
    positions_split2 = [positions_dict_split2[gene] for gene in common_genes]
    cor_spliced = compute_gene_correlation_between_splits(adata1.layers['spliced_original'][:,positions_split1],
                                                          adata2.layers['spliced_original'][:,positions_split2])
    cor_unspliced = compute_gene_correlation_between_splits(adata1.layers['unspliced_original'][:,positions_split1],
                                                            adata2.layers['unspliced_original'][:,positions_split2])
    Ngenes_spliced = len(cor_spliced[~np.isnan(cor_spliced)])
    Ngenes_unspliced = len(cor_spliced[~np.isnan(cor_spliced)])
    plt.clf()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
    # Plotting spliced
    axes[0].scatter(range(Ngenes_spliced), cor_spliced[~np.isnan(cor_spliced)],color='royalblue',alpha=0.4)
    axes[0].set_title("Correlation of gene expr between splits (spliced), N="+str(Ngenes_spliced))
    axes[0].set_xlabel("Genes")
    axes[0].set_ylabel("Correlation")
    # Plotting unspliced
    axes[1].scatter(range(Ngenes_unspliced), cor_unspliced[~np.isnan(cor_unspliced)],color='royalblue',alpha=0.4)
    axes[1].set_title("Correlation of gene expr between splits (unspliced), N="+str(Ngenes_unspliced))
    axes[1].set_xlabel("Genes")
    axes[1].set_ylabel("Correlation")
    # Adjusting layout to avoid overlap
    plt.tight_layout()
    plt.savefig(fig_folder+fig_path) 
    plt.clf()

######################################################
######################################################
## 3plots.py
######################################################
# velocity
### plot velocity on umap
def plot_velocities_scv_utv(adata_in,adata_raw,fig_folder,fig_info,dataset,method,color_label=None):
    if dataset=="ery":
        color_label = 'celltype'
    elif dataset=="pan":
        color_label = 'clusters'
    data_method = dataset+"_"+method
    # umapCompute
    scv.pl.velocity_embedding_stream(adata_in, basis='umap',color=color_label,save=fig_folder+"velocity/"+data_method+"_"+fig_info+"_umapCompute.png")
    # umapOriginal
    adata=adata_in.copy()
    adata.obsm['X_umap'] = adata_raw.obsm['X_umap'].copy()
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=color_label,save=fig_folder+"velocity/"+data_method+"_"+fig_info+"_umapOriginal.png")    

######################################################
# cosine similarity
### compute cosine similarity
def compute_cosine_similarity(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    if method=="scv":
        velo_genes_split1 = adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
        velo_genes_split2 = adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    common_genes_velocity = np.intersect1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Number of overlapped genes for velocity computation in splits = '+str(common_genes_velocity.shape[0])) 
    velo_df1 = pd.DataFrame(adata_split1.layers['velocity'], columns=adata_split1.var.index.tolist())
    velo_df2 = pd.DataFrame(adata_split2.layers['velocity'], columns=adata_split2.var.index.tolist())
    cos_sim = np.diag(cosine_similarity(velo_df1[common_genes_velocity],velo_df2[common_genes_velocity]))
    return cos_sim, common_genes_velocity.shape[0] # return cosine similarity and number of common genes in velocity computation

### plot cosine similarity
def plot_cosine_similarity(adata_split1,adata_split2,adata_total,adata_raw,dataset,method,fig_folder,text_x=None,text_y=None):
    cos_sim, Ngenes = compute_cosine_similarity(adata_split1,adata_split2,method)
    dataset_method = dataset+'_'+method
    # histogram
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='dimgray',color='powderblue') 
    max_frequency = np.max(counts)
    if text_x is None:
        text_x = np.quantile(cos_sim,[.05])[0]
    if text_y is None:
        text_y = max_frequency/2
    plt.axvline(np.mean(cos_sim), color='salmon', linestyle='dashed', linewidth=1.5) ## add mean
    plt.text(text_x,text_y,'mean='+str(np.round(np.mean(cos_sim),4))+', median='+str(np.round(np.median(cos_sim),4)), color='navy', fontsize=11)
    plt.xlabel('cosine similarity')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, '+dataset+'+'+method+', Ngenes='+str(Ngenes))
    plt.savefig(fig_folder+'cos_sim/'+dataset_method+'_cos_sim_hist.png')
    plt.clf()
    # umapCompute
    adata_total_plot = adata_total.copy()
    adata_total_plot.obs['cos_sim'] = cos_sim
    scv.pl.velocity_embedding_stream(adata_total_plot, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                     save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_umapCompute.png")
    scv.pl.scatter(adata_total_plot, color='cos_sim', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_scatter_umapCompute.png")
    # umapOriginal
    adata_total_plot = adata_total.copy()
    adata_total_plot.obs['cos_sim'] = cos_sim
    adata_total_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    scv.pl.velocity_embedding_stream(adata_total_plot, basis='umap',color="cos_sim",cmap='coolwarm',perc=[1, 100],
                                     save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_umapOriginal.png")
    scv.pl.scatter(adata_total_plot, color='cos_sim', cmap='coolwarm', perc=[1, 100],
                   save=fig_folder+"cos_sim/"+dataset_method+"_cos_sim_scatter_umapOriginal.png")
#### plot 2 by 1 figures, left: umap of cell development with velocity estimates, right: cosine similarity of velocities
def plot_cosine_similarity_withRef(adata_split1,adata_split2,adata_total,adata_raw,dataset,method,fig_folder,text_x=None,text_y=None):
    cos_sim, Ngenes = compute_cosine_similarity(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset+'_'+method
    celltype_label = None
    if dataset=="pan":
        celltype_label = "clusters"
    if dataset=="ery":
        celltype_label = "celltype"
    # umapCompute
    adata_plot = adata_total.copy()
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata_plot,color='cos_sim',cmap='coolwarm',perc=[1,100],ax=axs[1],legend_loc='none',
                   title='cosine similarity, '+dataset+'+'+method+', Ngenes='+str(Ngenes),frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+"cos_sim/"+dataset_method+"_cos_sim_withRef_umapCompute.png")
    plt.clf()
    # umapOriginal
    adata_plot = adata_total.copy()
    adata_plot.obsm['X_umap'] = adata_raw.obsm['X_umap']
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata_plot,color='cos_sim',cmap='coolwarm',perc=[1,100],ax=axs[1],legend_loc='none',
                   title='cosine similarity, '+dataset+'+'+method+', Ngenes='+str(Ngenes), frameon=False,size=100,alpha=0.3)
    plt.savefig(fig_folder+"cos_sim/"+dataset_method+"_cos_sim_withRef_umapOriginal.png")
    plt.clf()

######################################################
# plot pseudotime
def plot_pseudotime(adata_in,data_raw,fig_name,dataset,method):
    data_method = dataset+'_'+method
    fig_name = data_method+'_'+fig_name
    if not 'velocity_pseudotime' in adata_in.obs.columns: raise ValueError('No pseudotime information')
    celltype_label = "celltype"
    if dataset=="pan": celltype_label="clusters"
    # umapCompute
    adata = adata_in.copy()
    scv.pl.scatter(adata, color="velocity_pseudotime", color_map="gnuplot",title='velocity pseudotime, '+dataset+'+'+method,
                   save=fig_folder+'ptime/'+fig_name+'_ptime_umapCompute.png')
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color="velocity_pseudotime", color_map="gnuplot",title='velocity pseudotime, '+dataset+'+'+method)
    plt.savefig(fig_folder+"ptime/"+fig_name+'_ptime_withRef_umapCompute.png')
    plt.clf()
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = data_raw.obsm['X_umap'].copy()
    scv.pl.scatter(adata, color="velocity_pseudotime", color_map="gnuplot", title='velocity pseudotime, '+dataset+'+'+method,
                   save=fig_folder+'ptime/'+fig_name+'_ptime_umapOriginal.png')
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color="velocity_pseudotime", color_map="gnuplot",title='velocity pseudotime, '+dataset+'+'+method)
    plt.savefig(fig_folder+"ptime/"+fig_name+'_ptime_withRef_umapOriginal.png')
    plt.clf()

# plot pseudotime correlation
def ptime_correlation_scatter_plot(s1,s2,method,dataset,name,xlab,ylab,fig_folder):
    celltype_label = "celltype"
    if dataset == "pan": celltype_label = "clusters"
    cell_types = s1.obs[celltype_label]
    colors = dict(zip(s1.obs[celltype_label].cat.categories, s1.uns[celltype_label+'_colors']))
    #raw.obs['clusters'].cat.categories,raw.uns['clusters_colors']
    data_method = dataset+'_'+method
    df = pd.DataFrame({'split1':s1.obs['velocity_pseudotime'],'split2':s2.obs['velocity_pseudotime'],'cell_types':cell_types})
    corr = np.round(np.corrcoef(s1.obs['velocity_pseudotime'],s2.obs['velocity_pseudotime']),3)[0,1]
    print(corr)
    plt.figure(figsize=(7, 5))
    for category, color in colors.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors))
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime correlation '+name+', '+dataset+'+'+method+', corr='+str(corr)+')')
    plt.savefig(fig_folder+'ptime/'+data_method+"_pseudotimeCorr"+name+".png")
    plt.close()
