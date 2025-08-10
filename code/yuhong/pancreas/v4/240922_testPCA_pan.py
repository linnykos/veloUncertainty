import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import random

def create_df_with_gene_columns(adata,genes):
    df = pd.DataFrame(columns=genes)
    for gene in np.intersect1d(genes,adata.var.index):
        idx = np.where(adata.var.index==gene)[0][0]
        df[gene] = adata.layers['velocity'][:,idx]
    return df

## want to project df(1|2)_uni onto df(1|2)_intersect
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.metrics.pairwise import cosine_similarity

def divide_folds_by_celltype(adata,celltype_label,n_folds=10,seed=1530):
    celltypes = adata.obs[celltype_label].cat.categories
    folds = [[] for _ in range(n_folds)]
    random.seed(seed)
    np.random.seed(seed)
    for ct in celltypes:
        idx_ct = np.where(adata.obs[celltype_label]==ct)[0]
        fold_assn = np.array(random.choices(range(n_folds), k=len(idx_ct)))
        for i in range(n_folds):
            folds[i].extend(idx_ct[np.where(fold_assn == i)[0]].tolist())
    return folds

def cv_pca_by_celltype(adata,genes_intersect,celltype_label,components=[5],n_folds=10,seed=1530):
    folds = divide_folds_by_celltype(adata,celltype_label,n_folds,seed)
    v = adata.layers['velocity']
    mse_res = np.zeros((len(components), n_folds)) # each column for a fold
    R2_res = np.zeros((len(components), n_folds))
    for k in range(n_folds):
        print('fold '+str(k))
        idx_test = folds[k]
        idx_train = np.setdiff1d(np.arange(split1.n_obs), idx_test)
        train_uni = create_df_with_gene_columns(adata[idx_train,], list(set(adata.var.index[~np.isnan(v[0])])-set(genes_intersect)))
        test_uni = create_df_with_gene_columns(adata[idx_test,], list(set(adata.var.index[~np.isnan(v[0])])-set(genes_intersect)))
        train_intersect = create_df_with_gene_columns(adata[idx_train,], list(genes_intersect))
        test_intersect = create_df_with_gene_columns(adata[idx_test,], list(genes_intersect))
        for i in range(len(components)):
            mse, R2 = fit_pca_train_and_test(train_uni,test_uni,train_intersect,test_intersect,n_components=components[i])
            mse_res[i,k] = mse
            R2_res[i,k] = R2
    return mse_res, R2_res

# fit uni onto intersect
def fit_pca_train_and_test(train_uni,test_uni,train_intersect,test_intersect,n_components=10):
    pca = PCA(n_components)
    Z = pca.fit_transform(train_uni)   # (N_train, n_components)
    regression = LinearRegression()
    regression.fit(Z, train_intersect)  # Learn the mapping from PCA space to intersected gene space
    pred = regression.predict(pca.fit_transform(test_uni))
    mse = mean_squared_error(test_intersect, pred)  
    R2 = r2_score(test_intersect, pred)
    print('mse='+str(np.round(mse,5))+', R2='+str(np.round(R2,5)))
    return mse,R2

def fit_pca_all_data(v_uni,v_intersect,n_components=10):
    pca = PCA(n_components)
    Z = pca.fit_transform(v_uni)  
    regression = LinearRegression()
    regression.fit(Z, v_intersect)  # Learn the mapping from PCA space to intersected gene space
    pred = regression.predict(Z)
    return pred

def compute_cos_sim_with_pca(split1,split2,celltype_label=None,celltype=None,n_components=[6,6]):
    v1 = split1.layers['velocity']
    v2 = split2.layers['velocity']
    genes_intersect = np.intersect1d(split1.var.index[~np.isnan(v1[0])], split2.var.index[~np.isnan(v2[0])])
    df1_intersect = create_df_with_gene_columns(split1,list(genes_intersect))
    df2_intersect = create_df_with_gene_columns(split2,list(genes_intersect))
    df1_uni = create_df_with_gene_columns(split1,list(set(split1.var.index[~np.isnan(v1[0])])-set(genes_intersect)))
    df2_uni = create_df_with_gene_columns(split2,list(set(split2.var.index[~np.isnan(v2[0])])-set(genes_intersect)))
    if not celltype == None: 
        idx_ct = np.where(split1.obs[celltype_label]==celltype)[0]
    else:
        celltype = 'All'
        idx_ct = range(split1.n_obs)
    pred1 = fit_pca_all_data(v_uni=np.array(df1_uni)[idx_ct,],v_intersect=np.array(df1_intersect)[idx_ct,],
                         n_components=n_components[0])
    pred2 = fit_pca_all_data(v_uni=np.array(df2_uni)[idx_ct,],v_intersect=np.array(df2_intersect)[idx_ct,],
                         n_components=n_components[1])
    print(celltype+': cos_sim_uni after pca, cos_sim_intersect, cos_sim_add with predicted+intersection genes, (intersect+uni)/2')
    cos_sim_uni = np.diag(cosine_similarity(pred1, pred2))
    print(np.quantile(cos_sim_uni,[0.,.25,.5,.75,1.]))
    cos_sim_intersect = np.diag(cosine_similarity(np.array(df1_intersect)[idx_ct,], np.array(df2_intersect)[idx_ct,]))
    print(np.quantile(cos_sim_intersect,[0.,.25,.5,.75,1.]))
    cos_sim_add = np.diag(cosine_similarity(pred1+np.array(df1_intersect)[idx_ct,], pred2+np.array(df2_intersect)[idx_ct,]))
    print(np.quantile(cos_sim_add,[0.,.25,.5,.75,1.]))
    cos_sim_mean_of2 = (cos_sim_intersect + cos_sim_uni)/2
    print(np.quantile(cos_sim_mean_of2,[0.,.25,.5,.75,1.]))


dataset_long = 'pancreas'
dataset_short = 'pan'
method = 'utv'
outputAdded = '_outputAdded'
if method=='scv': outputAdded=''

split1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
split2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')


compute_cos_sim_with_pca(split1,split2,celltype_label='clusters',celltype='Ngn3 high EP')

compute_cos_sim_with_pca(split1,split2,celltype_label='clusters')


# ['Ductal', 'Ngn3 low EP', 'Ngn3 high EP', 'Pre-endocrine', 'Beta', 'Alpha', 'Delta', 'Epsilon']

def cv_pca_components_method(method,dataset_short):
    if dataset_short=='pan': 
        dataset_long = 'pancreas'
        celltype_label = 'clusters'
    elif dataset_short=='ery': 
        dataset_long = 'erythroid'
        celltype_label = 'celltype'
    outputAdded = '_outputAdded'
    if method=='scv': outputAdded=''
    split1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
    split2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')
    v1 = split1.layers['velocity']
    v2 = split2.layers['velocity']
    genes_intersect = np.intersect1d(split1.var.index[~np.isnan(v1[0])], split2.var.index[~np.isnan(v2[0])])
    mse_split1,R2_split1 = cv_pca_by_celltype(adata=split1,genes_intersect=genes_intersect,celltype_label=celltype_label,components=[i+1 for i in range(10)])
    mse_split2,R2_split2 = cv_pca_by_celltype(adata=split2,genes_intersect=genes_intersect,celltype_label=celltype_label,components=[i+1 for i in range(10)])
    mse_mean1 = np.mean(mse_split1,1)
    R2_mean1 = np.mean(R2_split1,1)
    mse_mean2 = np.mean(mse_split2,1)
    R2_mean2 = np.mean(R2_split2,1)
    return mse_mean1,R2_mean1,mse_mean2,R2_mean2

cv_pca_components_method(method='scv',dataset_short='pan')
    
dataset_long = 'pancreas'
dataset_short = 'pan'
method = 'velovi_woprep'
outputAdded = '_outputAdded'
if method=='scv': outputAdded=''

split1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
split2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')

v1 = split1.layers['velocity']
v2 = split2.layers['velocity']
genes_intersect = np.intersect1d(split1.var.index[~np.isnan(v1[0])], split2.var.index[~np.isnan(v2[0])])

mse_split1,R2_split1 = cv_pca_by_celltype(adata=split1,genes_intersect=genes_intersect,celltype_label='clusters',components=[i+1 for i in range(15)])
np.mean(mse_split1,1)
np.mean(R2_split1,1)

mse_split2,R2_split2 = cv_pca_by_celltype(adata=split2,genes_intersect=genes_intersect,celltype_label='clusters',components=[i+1 for i in range(15)])
np.mean(mse_split2,1)
np.mean(R2_split2,1)


compute_cos_sim_with_pca(split1,split2,celltype_label='celltype')

def compute_cos_sim_with_pca_method(method,dataset_short,n_components):
    if dataset_short=='pan': 
        dataset_long = 'pancreas'
        celltype_label = 'clusters'
    elif dataset_short=='ery': 
        dataset_long = 'erythroid'
        celltype_label = 'celltype'
    outputAdded = '_outputAdded'
    if method=='scv': outputAdded=''
    split1 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split1_v4'+outputAdded+'.h5ad')
    split2 = sc.read_h5ad('/Users/y2564li/Downloads/proj_scRNA/data/v4_'+dataset_long+'/seed317/'+method+'/adata_'+dataset_short+'_'+method+'_split2_v4'+outputAdded+'.h5ad')
    v1 = split1.layers['velocity']
    v2 = split2.layers['velocity']
    genes_intersect = np.intersect1d(split1.var.index[~np.isnan(v1[0])], split2.var.index[~np.isnan(v2[0])])
    compute_cos_sim_with_pca(split1,split2,celltype_label,celltype=None,n_components=n_components)

compute_cos_sim_with_pca_method(method='scv',dataset_short='ery',n_components=[6,6])
compute_cos_sim_with_pca_method(method='utv',dataset_short='ery',n_components=[6,8])
compute_cos_sim_with_pca_method(method='sct',dataset_short='ery',n_components=[7,9])
compute_cos_sim_with_pca_method(method='velovi',dataset_short='ery',n_components=[5,6])
compute_cos_sim_with_pca_method(method='velovi_woprep',dataset_short='ery',n_components=[5,7])


compute_cos_sim_with_pca_method(method='scv',dataset_short='pan',n_components=[7,6])
compute_cos_sim_with_pca_method(method='utv',dataset_short='pan',n_components=[11,13])
compute_cos_sim_with_pca_method(method='sct',dataset_short='pan',n_components=[7,9])
compute_cos_sim_with_pca_method(method='velovi',dataset_short='pan',n_components=[2,3])
compute_cos_sim_with_pca_method(method='velovi_woprep',dataset_short='pan',n_components=[2,1])

from scipy import stats
stats.spearmanr(split1.layers['latent_time_velovi'][:,0],split1.layers['latent_time_velovi'][:,1])