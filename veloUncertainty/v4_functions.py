import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
import scanpy as sc
import datetime

def read_data_v4(dataset_long,dataset_short,method,split_seed,data_version,allgenes=False,outputAdded=False):
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'
    if (allgenes==False):
        data_path = data_folder+'v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_'+data_version+'_v4'
    elif 'split' in data_version:
        data_path = data_folder+'v4_'+dataset_long+'/'+'seed'+str(split_seed)+'_'+dataset_short+'_split1_allgenes'
    elif data_version=='total':
        data_path = data_folder+'v4_'+dataset_long+'/'+dataset_short+'_total_allgenes'
    if outputAdded:
        data_path = data_path+'_outputAdded'
    print(data_path+'.h5ad')
    return sc.read_h5ad(data_path+'.h5ad')

def get_umap_sct(adata,umapOriginal=False,moments=True,velocity_graph=True):
    if moments==True:
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=30)
    if umapOriginal==True:
        adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    else:
        sc.tl.umap(adata)
    if velocity_graph==True: 
        scv.tl.velocity_graph(adata,n_jobs=8)

## print current time with a message
def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{message} at {current_time}")

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

def read_raw_adata(dataset):
    if 'ery' in dataset: 
        return sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad')
    elif ('pan' in dataset) and ('INC' in dataset):
        return sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Pancreas/endocrinogenesis_day15.h5ad')
    elif ('pan' in dataset) and (not 'INC' in dataset):
        return sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v2_pancreasINC/pancreasINC_total_allgenes.h5ad')
    elif dataset=='larry':
        return sc.read_h5ad('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_larry/larry.h5ad')

def plot_method_gene_corr(split1, split2, method, dataset, fig_folder, split_seed):
    if 'ery' in dataset: 
        dataset_short = 'ery'
        dataset_long = 'erythroid'
    elif ('pan' in dataset) and ('INC' in dataset):
        dataset_short = 'panINC'
        dataset_long = 'pancreasINC'
    elif ('pan' in dataset) and (not 'INC' in dataset):
        dataset_short = 'pan'
        dataset_long = 'pancreas'
    elif 'larry' in dataset:
        dataset_short = 'larry'
        dataset_long = 'larry'
    celltype_label = get_celltype_label(dataset_short)
    common_genes = np.intersect1d(np.array(split1.var.index[np.where(~np.isnan(split1.layers['velocity'][0]))]), np.array(split2.var.index[np.where(~np.isnan(split2.layers['velocity'][0]))]))
    gene_names_split1 = split1.var.index.copy()
    positions_dict_split1 = {gene: pos for pos, gene in enumerate(gene_names_split1)}
    positions_split1 = [positions_dict_split1[gene] for gene in common_genes]
    gene_names_split2 = split2.var.index.copy()
    positions_dict_split2 = {gene: pos for pos, gene in enumerate(gene_names_split2)}
    positions_split2 = [positions_dict_split2[gene] for gene in common_genes]
    cor_spliced = compute_gene_correlation_between_splits(split1.layers['spliced_original'][:,positions_split1],
                                                        split2.layers['spliced_original'][:,positions_split2])
    cor_unspliced = compute_gene_correlation_between_splits(split1.layers['unspliced_original'][:,positions_split1],
                                                            split2.layers['unspliced_original'][:,positions_split2])
    Ngenes_spliced = len(cor_spliced[~np.isnan(cor_spliced)])
    Ngenes_unspliced = len(cor_unspliced[~np.isnan(cor_unspliced)])
    print('Ngenes_spliced='+str(Ngenes_spliced)+', '+'Ngenes_unspliced='+str(Ngenes_unspliced)+', Ngenes_common_total='+str(len(common_genes)))
    overdisp_S = np.array(pd.read_csv('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'+dataset_short+'_overdisp_S.csv')['x'])
    overdisp_U = np.array(pd.read_csv('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'+dataset_short+'_overdisp_U.csv')['x'])
    df = pd.DataFrame(columns=['gene_names'])
    raw = read_raw_adata(dataset)
    df['gene_names'] = raw.var.index
    df['overdisps_S'] = overdisp_S
    df['overdisps_U'] = overdisp_U
    df = df[df['gene_names'].isin(common_genes)]
    df['gene_names'] = pd.Categorical(df['gene_names'], categories=common_genes, ordered=True)
    df = df.sort_values('gene_names')
    df['cor_spliced'] = cor_spliced
    df['cor_unspliced'] = cor_unspliced
    plt.clf()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))  # 1 row, 2 columns
    axes[0].scatter(df['overdisps_S'].clip(upper=50), df['cor_spliced'], color='royalblue',alpha=0.4)
    axes[0].set_title('Correlation btw splits,'+method+'+'+dataset_short+' (seed='+str(split_seed)+', spliced), N='+str(np.sum(~np.isnan(df['overdisps_S']+df['cor_spliced']))))
    axes[0].set_xlabel('overdispersion (clipped by 50)')
    axes[0].set_ylabel('Correlation')
    axes[1].scatter(df['overdisps_U'].clip(upper=50), df['cor_unspliced'],color='royalblue',alpha=0.4)
    axes[1].set_title('Correlation btw splits,'+method+'+'+dataset_short+' (seed='+str(split_seed)+', unspliced), N='+str(np.sum(~np.isnan(df['overdisps_U']+df['cor_unspliced']))))
    axes[1].set_xlabel('overdispersion (clipped by 50)')
    axes[1].set_ylabel('Correlation')
    plt.tight_layout()
    plt.savefig(fig_folder+method+'_'+dataset_short+'_gene_corr_overdisps.png') 


# compute cosine similarity
def compute_cosine_similarity_intersect(adata_split1,adata_split2,method):
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

def compute_cosine_similarity_union(adata_split1,adata_split2,method):
    velo_genes_split1 = adata_split1.var.index
    velo_genes_split2 = adata_split2.var.index
    velo_split1 = pd.DataFrame(adata_split1.layers['velocity'], columns=velo_genes_split1)
    velo_split2 = pd.DataFrame(adata_split2.layers['velocity'], columns=velo_genes_split2)
    if method=='scv':
        velo_genes_split1 = velo_genes_split1[~np.isnan(velo_split1.loc[0])] #adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
        velo_genes_split2 = velo_genes_split2[~np.isnan(velo_split2.loc[0])] #adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    union_genes_velo = np.union1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Size of the union of genes for velocity computation in splits = '+str(union_genes_velo.shape[0])) 
    Nrow = adata_split1.shape[0]
    velo_df1 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split1:
        velo_df1[gene] = velo_split1[gene]
    velo_df2 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split2:
        velo_df2[gene] = velo_split2[gene]
    cos_sim = np.diag(cosine_similarity(velo_df1,velo_df2))
    return cos_sim, union_genes_velo.shape[0]

################################################
# plot metrics
### helper
def get_celltype_label(dataset):
    celltype_label = None
    if 'ery' in dataset: celltype_label = 'celltype'
    elif 'pan' in dataset: celltype_label = 'clusters'
    elif 'larry' in dataset: celltype_label = 'state_info'
    return celltype_label

### helper
def get_basis_type(basis_type):
    if 'C' in basis_type or 'c' in basis_type: basis_type = 'Compute'
    elif 'O' in basis_type or 'o' in basis_type: basis_type = 'Original'
    return basis_type

### helper
def get_metric_color_and_title(metric):
    metric_color = ''
    metric_title = ''
    if 'cos' in metric: 
        metric_color = 'cos_sim'
        metric_title = 'Cosine similarity'
    elif 'conf' in metric: 
        metric_color = 'velocity_confidence'
        metric_title = 'Velocity confidence'
    return metric_color, metric_title

### helper
def plot_metric(adata,metric,plot_type,fig_folder,dataset,method,basis_type,split_seed,basis='umap',Ngenes=None,recompute=True,cmap='coolwarm'):
    metric_color, metric_title = get_metric_color_and_title(metric)
    basis_type = get_basis_type(basis_type)
    Ngenes_title = ''
    if not Ngenes==None: Ngenes_title = ', Ngenes='+str(Ngenes)
    if basis_type=='Original': scv.tl.velocity_graph(adata)
    if 'emb' in plot_type:
        scv.pl.velocity_embedding_stream(adata, basis=basis, color=metric_color, cmap=cmap, recompute=recompute,
                                         perc=[1, 100], title=metric_title+" "+dataset+'+'+method+Ngenes_title+', split_seed='+str(split_seed),
                                         save=fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_'+basis+basis_type+'.png')
    elif 'sca' in plot_type:
        scv.pl.scatter(adata, color=metric_color, cmap=cmap, title=metric_title+" "+dataset+'+'+method+Ngenes_title+', split_seed='+str(split_seed), perc=[0, 100], 
                       save=fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_'+basis+basis_type+'_scatter.png')


# plot cosine similarity
def plot_cosine_similarity(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,recompute=True,text_x=None,text_y=None):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    print('median cos_sim='+str(np.median(cos_sim))+', mean cos_sim='+str(np.mean(cos_sim)))
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset+'_'+method
    # histogram
    print('Plot histogram')
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(cos_sim, bins=30, edgecolor='gainsboro',color='powderblue') 
    max_frequency = np.max(counts)
    if text_x is None: text_x = np.quantile(cos_sim,[.05])[0]
    if text_y is None: text_y = max_frequency/5
    plt.axvline(np.mean(cos_sim), color='brown', linestyle='dashed', linewidth=1.5) ## add mean
    plt.axvline(np.median(cos_sim), color='peru', linestyle='dashed', linewidth=1.5) ## add median
    plt.text(text_x,text_y*2.5,'mean='+str(np.round(np.mean(cos_sim),4)), color='firebrick', fontsize=11)
    plt.text(text_x,text_y*3,'median='+str(np.round(np.median(cos_sim),4)), color='sienna', fontsize=11)
    plt.xlabel('cosine similarity')
    plt.ylabel('Frequency')
    plt.title('Histogram of cosine similarity, '+dataset+'+'+method+', Ngenes='+str(Ngenes)+', split_seed='+str(split_seed))
    plt.savefig(fig_folder+'metric/'+dataset_method+'_cos_sim_hist.png')
    plt.clf()
    # umapCompute
    print('Plot umapCompute')
    adata_total_plot = adata_total.copy()
    plot_metric(adata=adata_total_plot,metric='cos',plot_type='emb',fig_folder=fig_folder,dataset=dataset,method=method,Ngenes=Ngenes,split_seed=split_seed,basis_type='Compute',basis='umap',recompute=recompute)
    plot_metric(adata=adata_total_plot,metric='cos',plot_type='scat',fig_folder=fig_folder,dataset=dataset,method=method,Ngenes=Ngenes,split_seed=split_seed,basis_type='Compute',basis='umap',recompute=recompute)
    # umapOriginal
    adata_total_plot = adata_total.copy()
    adata_total_plot.obsm['X_umap'] = adata_total_plot.obsm['X_umapOriginal']
    scv.tl.velocity_graph(adata_total_plot)
    print('Plot umapOriginal')
    plot_metric(adata=adata_total_plot,metric='cos',plot_type='emb',fig_folder=fig_folder,dataset=dataset,method=method,basis_type='Orig',basis='umap',Ngenes=Ngenes,split_seed=split_seed)
    plot_metric(adata=adata_total_plot,metric='cos',plot_type='scat',fig_folder=fig_folder,dataset=dataset,method=method,basis_type='Orig',basis='umap',Ngenes=Ngenes,split_seed=split_seed)

### helper
def plot_metric_withRef(adata,metric,dataset,method,fig_folder,basis_type,split_seed,celltype_label=None,Ngenes=None,recompute=True,basis='umap'):
    metric_color, metric_title = get_metric_color_and_title(metric)
    basis_type = get_basis_type(basis_type)
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    Ngenes_title = ''
    if not Ngenes==None: Ngenes_title = ', Ngenes='+str(Ngenes)
    if basis_type=='Original':
        scv.tl.velocity_graph(adata,n_jobs=8)
    fig,axs = plt.subplots(ncols=2, nrows=1, figsize=(11,4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata,basis=basis,color=celltype_label,ax=axs[0],legend_loc='on data',recompute=recompute,
                                     title="Velocity "+dataset+'+'+method, frameon=False, size=100, alpha=0.5)
    scv.pl.scatter(adata,color=metric_color,cmap='coolwarm',perc=[0,100],ax=axs[1],legend_loc='none',
                   title=metric_title+" "+dataset+'+'+method+Ngenes_title+', split_seed='+str(split_seed),frameon=False,size=100,alpha=0.2)
    plt.savefig(fig_folder+"metric/"+dataset+'_'+method+'_'+metric_color+'_withRef_'+basis+basis_type+'.png')
    plt.clf()

#### plot 2 by 1 figures, left: umap of cell development with velocity estimates, right: cosine similarity of velocities
def plot_cosine_similarity_withRef(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,recompute=True,celltype_label=None):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    # umapCompute
    adata_plot = adata_total.copy()
    plot_metric_withRef(adata=adata_plot,metric='cos',dataset=dataset,method=method,fig_folder=fig_folder,basis_type='C',split_seed=split_seed,celltype_label=celltype_label,Ngenes=Ngenes,recompute=recompute,basis='umap')
    # umapOriginal
    adata_plot = adata_total.copy()
    adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umapOriginal']
    scv.tl.velocity_graph(adata_plot,n_jobs=8)
    print('Plot umapOriginal')
    plot_metric_withRef(adata=adata_plot,metric='cos',dataset=dataset,method=method,fig_folder=fig_folder,basis_type='Orig',split_seed=split_seed,celltype_label=celltype_label,Ngenes=Ngenes,recompute=recompute,basis='umap')
    

######################################################
## helper: plot velo_conf
def plot_veloConf_and_cosSim_helper(adata_total,dataset,method,fig_folder,umapOriginal,Ngenes,split_seed,recompute=True,celltype_label=None):
    adata_plot = adata_total.copy()
    if celltype_label==None:  celltype_label = get_celltype_label(dataset)
    data_method = dataset+'_'+method
    fig_umap = "umapCompute"
    if umapOriginal==True:
        adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umapOriginal']
        scv.tl.velocity_graph(adata_plot,n_jobs=8)
        fig_umap = "umapOriginal"
    vmin = np.min([0, np.min(adata_plot.obs['cos_sim'])-1e-5, np.min(adata_plot.obs['velocity_confidence'])-1e-5])
    vmax = np.max([np.max(adata_plot.obs['cos_sim'])+1e-5, np.max(adata_plot.obs['velocity_confidence'])+1e-5, 1])
    Ngenes_conf = np.sum(~np.isnan(adata_total.layers['velocity'][0]))
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 4))  # figsize=(horizontal, vertical)
    scv.pl.velocity_embedding_stream(adata_plot, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     recompute=recompute,frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata_plot,c='velocity_confidence',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[1],legend_loc='none',
                   title='Velocity confidence, '+dataset+'+'+method+', Ngenes='+str(Ngenes_conf),frameon=False,size=100,alpha=0.3)
    scv.pl.scatter(adata_plot,color='cos_sim',cmap='coolwarm',vmin=vmin,vmax=vmax,ax=axs[2],legend_loc='none',frameon=False,
                   title='Velocity cosine similarity, '+dataset+'+'+method+', Ngenes='+str(Ngenes)+', split_seed='+str(split_seed),size=100,alpha=0.3)
    plt.savefig(fig_folder+"metric/"+data_method+"_veloConf_and_cosSim_"+fig_umap+".png")
    plt.clf()
    plot_metric(adata_plot,metric='conf',plot_type='scat',fig_folder=fig_folder,dataset=dataset,method=method,basis_type=fig_umap,split_seed=split_seed,Ngenes=Ngenes)
  
def plot_veloConf_and_cosSim(adata_total,adata_split1,adata_split2,dataset,method,fig_folder,split_seed,recompute=True):
    adata_plot = adata_total.copy()
    if (not 'velocity_confidence' in adata_plot.obs.columns):
        scv.tl.velocity_confidence(adata_plot)
    cos_sim,Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method=method)
    adata_plot.obs['cos_sim'] = cos_sim
    # umapCompute
    print('Plot umapCompute')
    plot_veloConf_and_cosSim_helper(adata_plot,dataset,method,fig_folder,umapOriginal=False,Ngenes=Ngenes,split_seed=split_seed,recompute=recompute)
    # umapOriginal
    print('Plot umapOriginal')
    plot_veloConf_and_cosSim_helper(adata_plot,dataset,method,fig_folder,umapOriginal=True,Ngenes=Ngenes,split_seed=split_seed,recompute=recompute)


def plot_veloConf_hist(adata_total,dataset,method,fig_folder,split_seed,text_x=None,text_y=None):
    adata_plot = adata_total.copy()
    if not 'velocity_confidence' in adata_plot.obs.columns:
        scv.tl.velocity_confidence(adata_plot)
    velo_conf = adata_plot.obs['velocity_confidence']
    Ngenes = len(adata_plot.layers['velocity'][0]) - np.sum(np.isnan(adata_plot.layers['velocity'][0]))
    # histogram
    plt.clf()
    plt.figure(figsize=(7, 5))
    counts, bins, patches = plt.hist(velo_conf, bins=30, edgecolor='gainsboro',color='powderblue') 
    max_frequency = np.max(counts)
    if text_x is None: text_x = np.quantile(velo_conf,[.05])[0]
    if text_y is None: text_y = max_frequency/5
    plt.axvline(np.mean(velo_conf), color='brown', linestyle='dashed', linewidth=1.5) ## add mean
    plt.axvline(np.median(velo_conf), color='peru', linestyle='dashed', linewidth=1.5) ## add median
    plt.text(text_x,text_y*2.5,'mean='+str(np.round(np.mean(velo_conf),4)),color='firebrick',fontsize=11)
    plt.text(text_x,text_y*3,'median='+str(np.round(np.median(velo_conf),4)),color='sienna',fontsize=11)
    plt.xlabel('Velocity confidence')
    plt.ylabel('Frequency')
    plt.title('Histogram of velocity confidence, '+dataset+'+'+method+', Ngenes='+str(Ngenes)+', split_seed='+str(split_seed))
    plt.savefig(fig_folder+'metric/'+dataset+'_'+method+'_veloConf_hist.png')
    plt.clf()



######################################################
# plot pseudotime
def plot_pseudotime(adata_in,data_version,dataset,method,fig_folder,split_seed,recompute=True,celltype_label=None,ptime_label='velocity_pseudotime'):
    fig_title = data_version+', split_seed='+str(split_seed)
    if not ptime_label in adata_in.obs.columns: raise ValueError('No pseudotime information')
    if celltype_label == None: celltype_label = get_celltype_label(dataset)
    # umapCompute
    adata = adata_in.copy()
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     recompute=recompute,frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color=ptime_label, color_map="gnuplot",title='pseudotime, '+dataset+'+'+method+' '+fig_title)
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+'_ptime_withRef_'+data_version+'_umapCompute.png')
    plt.clf()
    plt.clf()
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    scv.tl.velocity_graph(adata)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     recompute=recompute,frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color=ptime_label, color_map="gnuplot",title='pseudotime, '+dataset+'+'+method+' '+fig_title)
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+'_ptime_withRef_'+data_version+'_umapOriginal.png')
    plt.clf()


# plot pseudotime correlation
def ptime_correlation_scatter_spearman(s1,s2,method,dataset,name,xlab,ylab,fig_folder,time_label,split_seed,celltype_label=None,alpha=.3):
    from scipy.stats import spearmanr
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    cell_types = s1.obs[celltype_label]
    colors = dict(zip(s1.obs[celltype_label].cat.categories, s1.uns[celltype_label+'_colors']))
    time_type = 'pseudotime'
    df = pd.DataFrame({'split1':s1.obs[time_label],'split2':s2.obs[time_label],'cell_types':cell_types})
    corr = np.round(spearmanr(s1.obs[time_label], s2.obs[time_label]).correlation,3)
    print(corr)
    plt.figure(figsize=(7, 5))
    for category, color in colors.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors),alpha=alpha)
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime correlation '+name+', '+dataset+'+'+method+', split_seed='+str(split_seed)+' (corr='+str(corr)+')')
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+'_ptime_SpearmanCorr_'+name+'.png')
    plt.close()

# plot latent time
def plot_latent_time(adata_in,data_version,dataset,method,fig_folder,split_seed,recompute=True,celltype_label=None,time_label='latent_time'):
    data_method = dataset+'_'+method
    fig_title = data_version+', split_seed='+str(split_seed)
    data_version = data_method+'_'+data_version
    if not time_label in adata_in.obs.columns: raise ValueError('No latent_time information')
    if celltype_label == None: celltype_label = get_celltype_label(dataset)
    # umapCompute
    adata = adata_in.copy()
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     recompute=recompute,frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color=time_label, color_map="gnuplot",title='Latent time, '+dataset+'+'+method+' '+fig_title)
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+"_lattime_withRef_"+data_version+'_umapCompute.png')
    plt.clf()
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    scv.tl.velocity_graph(adata)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     recompute=recompute,frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color=time_label, color_map="gnuplot",title='Latent time, '+dataset+'+'+method+' '+fig_title)
    plt.savefig(fig_folder+"ptime/"+dataset+'_'+method+"_lattime_withRef_"+data_version+'_umapOriginal.png')
    plt.clf()

def latent_time_correlation_scatter_spearman(s1,s2,method,dataset,name,xlab,ylab,fig_folder,split_seed,celltype_label=None,alpha=.3):
    from scipy.stats import spearmanr
    time_label='latent_time'
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    cell_types = s1.obs[celltype_label]
    colors = dict(zip(s1.obs[celltype_label].cat.categories, s1.uns[celltype_label+'_colors']))
    df = pd.DataFrame({'split1':s1.obs[time_label],'split2':s2.obs[time_label],'cell_types':cell_types})
    corr = np.round(spearmanr(s1.obs[time_label], s2.obs[time_label]).correlation,3)
    print(corr)
    plt.figure(figsize=(7, 5))
    for category, color in colors.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors),alpha=alpha)
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Latent time correlation '+name+', '+dataset+'+'+method+', split_seed='+str(split_seed)+' (corr='+str(corr)+')')
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+"_latentSpearmanCorr_"+name+".png")
    plt.close()



################
## cosine similarity by celltype
def plot_cosine_similarity_hist_by_celltype_pancreas(adata_split1,adata_split2,adata_total,dataset,method,fig_folder):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset+'_'+method
    celltype_label = 'clusters' # 8 labels
    plt.clf()
    cell_cat = adata_total.obs[celltype_label].cat.categories
    fig,axs = plt.subplots(ncols=4, nrows=2, figsize=(25,10))
    axs = axs.ravel()
    for idx,ax in enumerate(axs):
        celltype = cell_cat[idx]
        cos_sim_celltype = cos_sim[adata_total.obs[celltype_label].array==celltype]
        Ncells = cos_sim_celltype.shape[0]
        counts, bins, patches = ax.hist(cos_sim_celltype, bins=20, edgecolor='gainsboro',color='powderblue') 
        max_frequency = np.max(counts)
        ax.axvline(np.mean(cos_sim_celltype), color='brown', linestyle='dashed', linewidth=1.5)
        ax.axvline(np.median(cos_sim_celltype), color='peru', linestyle='dashed', linewidth=1.5)
        text_x = np.quantile(cos_sim_celltype,[.0])[0]
        text_y = max_frequency/5
        ax.text(text_x,text_y*3,'mean='+str(np.round(np.mean(cos_sim_celltype),4)), color='firebrick', fontsize=11)
        ax.text(text_x,text_y*2,'median='+str(np.round(np.median(cos_sim_celltype),4)), color='sienna', fontsize=11)
        ax.set_xlabel('cosine similarity')
        ax.set_ylabel('Frequency')
        ax.set_title(celltype+' (Ncells='+str(Ncells)+'), '+dataset+'+'+method)
    plt.savefig(fig_folder+'metric/'+dataset_method+'_cos_sim_hist_byCelltype.png')
    plt.clf()

def plot_cosine_similarity_hist_by_celltype_pancreasINC(adata_split1,adata_split2,adata_total,dataset,method,fig_folder):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset+'_'+method
    celltype_label = 'clusters' # 7 labels
    plt.clf()
    cell_cat = adata_total.obs[celltype_label].cat.categories
    fig,axs = plt.subplots(ncols=4, nrows=2, figsize=(25,10))
    axs = axs.ravel()
    for idx,ax in enumerate(axs):
        if idx==7: break
        celltype = cell_cat[idx]
        cos_sim_celltype = cos_sim[adata_total.obs[celltype_label].array==celltype]
        Ncells = cos_sim_celltype.shape[0]
        counts, bins, patches = ax.hist(cos_sim_celltype, bins=20, edgecolor='gainsboro',color='powderblue') 
        max_frequency = np.max(counts)
        ax.axvline(np.mean(cos_sim_celltype), color='brown', linestyle='dashed', linewidth=1.5)
        ax.axvline(np.median(cos_sim_celltype), color='peru', linestyle='dashed', linewidth=1.5)
        text_x = np.quantile(cos_sim_celltype,[.0])[0]
        text_y = max_frequency/5
        ax.text(text_x,text_y*3,'mean='+str(np.round(np.mean(cos_sim_celltype),4)), color='firebrick', fontsize=11)
        ax.text(text_x,text_y*2,'median='+str(np.round(np.median(cos_sim_celltype),4)), color='sienna', fontsize=11)
        ax.set_xlabel('cosine similarity')
        ax.set_ylabel('Frequency')
        ax.set_title(celltype+' (Ncells='+str(Ncells)+'), '+dataset+'+'+method)
    plt.savefig(fig_folder+'metric/'+dataset_method+'_cos_sim_hist_byCelltype.png')
    plt.clf()

"""
def plot_cosine_similarity_hist_by_celltype(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,celltype_label_larry='state_info'):
    if dataset=='ery':
        return plot_cosine_similarity_hist_by_celltype_erythroid(adata_split1,adata_split2,adata_total,dataset,method,fig_folder)
    elif dataset=='pan':
        return plot_cosine_similarity_hist_by_celltype_pancreas(adata_split1,adata_split2,adata_total,dataset,method,fig_folder)
    elif dataset=='panINC':
        return plot_cosine_similarity_hist_by_celltype_pancreasINC(adata_split1,adata_split2,adata_total,dataset,method,fig_folder)
    elif dataset=='larry':
        return plot_cosine_similarity_hist_by_celltype_larry(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,celltype_label=celltype_label_larry)
"""

def plot_cosine_similarity_hist_by_celltype(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,celltype_label=None):
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset+'_'+method
    #celltype_label = 'state_info' # 11 labels
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    plt.clf()
    cell_cat = adata_total.obs[celltype_label].cat.categories
    nrows = len(cell_cat)%4
    fig,axs = plt.subplots(ncols=4, nrows=nrows, figsize=(25,nrows*5))
    axs = axs.ravel()
    for idx,ax in enumerate(axs):
        if idx==len(cell_cat): break
        celltype = cell_cat[idx]
        cos_sim_celltype = cos_sim[adata_total.obs[celltype_label].array==celltype]
        Ncells = cos_sim_celltype.shape[0]
        counts, bins, patches = ax.hist(cos_sim_celltype, bins=20, edgecolor='gainsboro',color='powderblue') 
        max_frequency = np.max(counts)
        ax.axvline(np.mean(cos_sim_celltype), color='brown', linestyle='dashed', linewidth=1.5)
        ax.axvline(np.median(cos_sim_celltype), color='peru', linestyle='dashed', linewidth=1.5)
        text_x = np.quantile(cos_sim_celltype,[.0])[0]
        text_y = max_frequency/5
        ax.text(text_x,text_y*3,'mean='+str(np.round(np.mean(cos_sim_celltype),4)), color='firebrick', fontsize=11)
        ax.text(text_x,text_y*2,'median='+str(np.round(np.median(cos_sim_celltype),4)), color='sienna', fontsize=11)
        ax.set_xlabel('cosine similarity')
        ax.set_ylabel('Frequency')
        ax.set_title(celltype+' (Ncells='+str(Ncells)+'), '+dataset+'+'+method)
    plt.savefig(fig_folder+'metric/'+dataset_method+'_cos_sim_byCelltype_hist.png')
    plt.clf()

def plot_cosine_similarity_hist_by_celltype(adata_split1,adata_split2,adata_total,dataset,method,fig_folder,split_seed,celltype_label=None):
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1,adata_split2,method)
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset+'_'+method
    plt.clf()
    cell_cat = adata_total.obs[celltype_label].cat.categories
    fig,axs = plt.subplots(ncols=4, nrows=3, figsize=(25,16))
    axs = axs.ravel()
    for idx,ax in enumerate(axs):
        if idx==len(cell_cat): break
        celltype = cell_cat[idx]
        cos_sim_celltype = cos_sim[adata_total.obs[celltype_label].array==celltype]
        Ncells = cos_sim_celltype.shape[0]
        counts, bins, patches = ax.hist(cos_sim_celltype, bins=20, edgecolor='gainsboro',color='powderblue') 
        max_frequency = np.max(counts)
        ax.axvline(np.mean(cos_sim_celltype), color='brown', linestyle='dashed', linewidth=1.5)
        ax.axvline(np.median(cos_sim_celltype), color='peru', linestyle='dashed', linewidth=1.5)
        text_x = np.quantile(cos_sim_celltype,[.0])[0]
        text_y = max_frequency/5
        ax.text(text_x,text_y*3,'mean='+str(np.round(np.mean(cos_sim_celltype),4)), color='firebrick', fontsize=11)
        ax.text(text_x,text_y*2,'median='+str(np.round(np.median(cos_sim_celltype),4)), color='sienna', fontsize=11)
        ax.set_xlabel('cosine similarity')
        ax.set_ylabel('Frequency')
        ax.set_title(celltype+' (Ncells='+str(Ncells)+'), '+dataset+'+'+method+', split_seed='+str(split_seed))
    plt.savefig(fig_folder+'metric/'+dataset_method+'_cos_sim_byCelltype_hist.png')
    plt.clf()

def plot_cosine_similarity_boxplot_by_celltype(adata_split1, adata_split2, adata_total, dataset, method, fig_folder,split_seed, celltype_label=None):
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    cos_sim, Ngenes = compute_cosine_similarity_union(adata_split1, adata_split2, method)
    adata_total.obs['cos_sim'] = cos_sim
    dataset_method = dataset + '_' + method
    plt.clf()
    cell_cat = adata_total.obs[celltype_label].cat.categories
    fig, ax = plt.subplots(figsize=(len(cell_cat)*1.6, 9))
    data_to_plot = [cos_sim[adata_total.obs[celltype_label].array==celltype] for celltype in cell_cat]
    # Create boxplot
    ax.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'),
               medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'),
               capprops=dict(color='rosybrown'), 
               flierprops=dict(marker='o', color='rosybrown', alpha=0.5, markerfacecolor='aliceblue', markeredgecolor='rosybrown'))
    counts = [len(data) for data in data_to_plot]
    x_labels = [f'{cell_cat[i]} (n={counts[i]})' for i in range(len(cell_cat))]
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=12)
    ax.set_xlabel('Cell Types')
    ax.set_ylabel('Cosine Similarity')
    ax.set_ylim(-1,1)
    ax.set_title(f'Cosine Similarity by Cell Type ({dataset} + {method})'+', Ngenes='+str(Ngenes)+', split_seed='+str(split_seed))
    plt.tight_layout()
    plt.savefig(fig_folder+'metric/'+dataset_method+'_cos_sim_byCelltype_boxplot.png')
    plt.clf()


def plot_velo_conf_boxplot_by_celltype(adata_plot,dataset,method,fig_folder,split_seed,celltype_label=None):
    if celltype_label==None: celltype_label = get_celltype_label(dataset)
    if (not method=='sct') and (not 'velocity_confidence' in adata_plot.obs.columns):
        scv.tl.velocity_confidence(adata_plot)
    velo_conf = adata_plot.obs['velocity_confidence']
    Ngenes = len(adata_plot.layers['velocity'][0]) - np.sum(np.isnan(adata_plot.layers['velocity'][0]))
    dataset_method = dataset + '_' + method
    plt.clf()
    cell_cat = adata_plot.obs[celltype_label].cat.categories
    fig, ax = plt.subplots(figsize=(len(cell_cat)*1.6, 9))
    data_to_plot = [velo_conf[adata_plot.obs[celltype_label].array==celltype] for celltype in cell_cat]
    # Create boxplot
    ax.boxplot(data_to_plot, patch_artist=True, boxprops=dict(facecolor='lightsteelblue', color='rosybrown'),
               medianprops=dict(color='rosybrown'), whiskerprops=dict(color='rosybrown'), capprops=dict(color='rosybrown'), 
               flierprops=dict(marker='o', color='rosybrown', alpha=0.5, markerfacecolor='aliceblue', markeredgecolor='rosybrown'))
    counts = [len(data) for data in data_to_plot]
    x_labels = [f'{cell_cat[i]} (n={counts[i]})' for i in range(len(cell_cat))]
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=12)
    ax.set_xlabel('Cell Types')
    ax.set_ylabel('Velocity Confidence')
    ax.set_ylim(-1,1)
    ax.set_title(f'Velocity Confidence by Cell Type ({dataset} + {method})'+', Ngenes='+str(Ngenes)+', split_seed='+str(split_seed))
    plt.tight_layout()
    plt.savefig(fig_folder+'metric/'+dataset_method+'_veloConf_byCelltype_boxplot.png')
    plt.clf()



def compute_cosine_similarity_shuffled(split1,split2,method,seed=1514,Niter=100):
    from sklearn.metrics.pairwise import cosine_similarity
    velo_genes_split1 = split1.var.index
    velo_genes_split2 = split2.var.index
    velo_split1 = pd.DataFrame(split1.layers['velocity'], columns=velo_genes_split1)
    velo_split2 = pd.DataFrame(split2.layers['velocity'], columns=velo_genes_split2)
    if method=='scv':
        velo_genes_split1 = velo_genes_split1[~np.isnan(velo_split1.loc[0])] #adata_split1.var.index[~np.isnan(adata_split1.layers['velocity'][0])]
        velo_genes_split2 = velo_genes_split2[~np.isnan(velo_split2.loc[0])] #adata_split2.var.index[~np.isnan(adata_split2.layers['velocity'][0])]
    union_genes_velo = np.union1d(np.array(velo_genes_split1), np.array(velo_genes_split2))
    print('Size of the union of genes for velocity computation in splits = '+str(union_genes_velo.shape[0])) 
    Nrow = split1.shape[0]
    velo_df1 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split1: velo_df1[gene] = velo_split1[gene]
    velo_df2 = pd.DataFrame(0, index=range(Nrow), columns=union_genes_velo)
    for gene in velo_genes_split2: velo_df2[gene] = velo_split2[gene]
    cos_sim = np.diag(cosine_similarity(velo_df1, velo_df2))
    print(np.round(np.mean(cos_sim),4), np.round(np.median(cos_sim),4))
    v2s_mean = []
    v2s_median = []
    np.random.seed(seed)
    for i in range(Niter+1):
        if i % 10==0: print(i)
        v2s = velo_df2.sample(frac=1)
        cos_sim_s = np.diag(cosine_similarity(velo_df1, v2s))
        v2s_mean.append(np.round(np.mean(cos_sim_s),4))
        v2s_median.append(np.round(np.median(cos_sim_s),4))
    return v2s_mean,v2s_median

def plot_pseudotime_diffusion(adata_in,data_version,dataset,method,fig_folder,split_seed,recompute=True,celltype_label=None,ptime_label='velocity_pseudotime'):
    fig_title = data_version+', split_seed='+str(split_seed)
    if not ptime_label in adata_in.obs.columns: raise ValueError('No pseudotime information')
    if celltype_label == None: celltype_label = get_celltype_label(dataset)
    # umapCompute
    adata = adata_in.copy()
    scv.tl.velocity_pseudotime(adata,use_velocity_graph=False)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     recompute=recompute,frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color=ptime_label, color_map="gnuplot",title='pseudotime, '+dataset+'+'+method+' '+fig_title)
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+'_ptime_withRef_diffusion_'+data_version+'_umapCompute.png')
    plt.clf()
    plt.clf()
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata,use_velocity_graph=False)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 5))
    scv.pl.velocity_embedding_stream(adata, basis='umap',color=celltype_label,ax=axs[0],legend_loc='on data',
                                     recompute=recompute,frameon=False,size=100,alpha=0.5)
    scv.pl.scatter(adata,ax=axs[1], color=ptime_label, color_map="gnuplot",title='pseudotime, '+dataset+'+'+method+' '+fig_title)
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+'_ptime_withRef_diffusion_'+data_version+'_umapOriginal.png')
    plt.clf()

def ptime_diffusion_correlation_scatter_spearman(s1,s2,method,dataset,name,xlab,ylab,fig_folder,time_label,split_seed,celltype_label=None,alpha=.3):
    scv.tl.velocity_pseudotime(s1,use_velocity_graph=False)
    scv.tl.velocity_pseudotime(s2,use_velocity_graph=False)
    from scipy.stats import spearmanr
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    cell_types = s1.obs[celltype_label]
    colors = dict(zip(s1.obs[celltype_label].cat.categories, s1.uns[celltype_label+'_colors']))
    time_type = 'pseudotime'
    df = pd.DataFrame({'split1':s1.obs[time_label],'split2':s2.obs[time_label],'cell_types':cell_types})
    corr = np.round(spearmanr(s1.obs[time_label], s2.obs[time_label]).correlation,3)
    print(corr)
    plt.figure(figsize=(7, 5))
    for category, color in colors.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors),alpha=alpha)
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime correlation '+name+', '+dataset+'+'+method+', split_seed='+str(split_seed)+' (corr='+str(corr)+')')
    plt.savefig(fig_folder+'ptime/'+dataset+'_'+method+'_ptime_SpearmanCorr_diffusion_'+name+'.png')
    plt.close()

def print_ptime_corr_by_celltype(split1,split2,total,dataset,ptime_label='velocity_pseudotime'):
    celltype_label = get_celltype_label(dataset)
    for i in range(len(total.obs[celltype_label].cat.categories)):
        ct = total.obs[celltype_label].cat.categories[i]
        print(ct, np.round(np.corrcoef(split1[split1.obs[celltype_label]==ct].obs[ptime_label],
                            split2[split2.obs[celltype_label]==ct].obs[ptime_label])[0,1],4))
