import scvelo as scv
import scanpy as sc
import datetime

import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import get_celltype_label

def print_message_with_time(message):
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'{message} at {current_time}')

def scv_compute_velocity(adata,dataset):
    if ('ery' in dataset): scv_compute_velocity_ery(adata)
    else:
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)
        sc.tl.umap(adata)
        scv.tl.recover_dynamics(adata,n_jobs=8)
        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata)

def scv_compute_velocity_ery(adata):
    import bbknn
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    ### batch correction
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    print("Batch correction done!")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40) 
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata,n_jobs=8)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

def plot_velocity_scv_utv(adata_in,fig_folder,fig_info,dataset,method,split_seed,recompute=True,celltype_label=None,basis='umap'):
    if celltype_label==None: celltype_label=get_celltype_label(dataset)
    data_method = dataset+"_"+method
    # umapCompute
    scv.pl.velocity_embedding_stream(adata_in, basis=basis,color=celltype_label,recompute=recompute,
                                     title='Velocity '+dataset+'+'+method+' '+fig_info+' (split_seed='+str(split_seed)+')',
                                     save=fig_folder+"velocity/"+data_method+"_"+fig_info+'_'+basis+"Compute.png")
    # umapOriginal
    adata = adata_in.copy()
    adata.obsm['X_umap'] = adata.obsm['X_umapOriginal'].copy()
    scv.pl.velocity_embedding_stream(adata, basis=basis,color=celltype_label,recompute=recompute,
                                     title='Velocity '+dataset+'+'+method+' '+fig_info+' (split_seed='+str(split_seed)+')',
                                     save=fig_folder+"velocity/"+data_method+"_"+fig_info+'_'+basis+"Original.png")    

