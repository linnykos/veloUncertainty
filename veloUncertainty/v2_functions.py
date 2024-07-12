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

## 3plots.py
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
