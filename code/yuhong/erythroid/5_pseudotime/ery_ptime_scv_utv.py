import scvelo as scv
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bbknn


def run_scv(adata):
    scv.pp.normalize_per_cell(adata)
    scv.pp.log1p(adata)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    ## batch correction
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    adata.X = adata.X.toarray()
    bbknn.ridge_regression(adata, batch_key='sample', confounder_key='celltype')
    sc.tl.pca(adata)
    bbknn.bbknn(adata, batch_key='sequencing.batch')
    ## batch correction done
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    return adata

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/"

adata_total = scv.read(data_folder+'erythroid_seed317_total_seurat.h5ad')
## total
run_scv(adata_total)
adata_total.write(filename=data_folder+'erythroid_seed317_total.h5ad')

## seed=317
adata_split1_seed317 = scv.read(data_folder+'erythroid_seed317_split1_seurat.h5ad')
adata_split2_seed317 = scv.read(data_folder+'erythroid_seed317_split2_seurat.h5ad')
run_scv(adata_split1_seed317)
run_scv(adata_split2_seed317)
adata_split1_seed317.write(filename=data_folder+'erythroid_seed317_split1.h5ad')
adata_split2_seed317.write(filename=data_folder+'erythroid_seed317_split2.h5ad')

np.corrcoef(adata_split1_seed317.obs['velocity_pseudotime'],adata_split2_seed317.obs['velocity_pseudotime'])
## -0.4035575
np.corrcoef(adata_split1_seed317.obs['velocity_pseudotime'],adata_total.obs['velocity_pseudotime']) # 0.63599917
np.corrcoef(adata_total.obs['velocity_pseudotime'],adata_split2_seed317.obs['velocity_pseudotime']) # 0.12670515


## seed=320
adata_split1_seed320 = scv.read(data_folder+'erythroid_seed320_split1_seurat.h5ad')
adata_split2_seed320 = scv.read(data_folder+'erythroid_seed320_split2_seurat.h5ad')
run_scv(adata_split1_seed320)
run_scv(adata_split2_seed320)
adata_split1_seed320.write(filename=data_folder+'erythroid_seed320_split1.h5ad')
adata_split2_seed320.write(filename=data_folder+'erythroid_seed320_split2.h5ad')

np.corrcoef(adata_split1_seed320.obs['velocity_pseudotime'],adata_split2_seed320.obs['velocity_pseudotime'])
## 0.04919378
np.corrcoef(adata_split1_seed320.obs['velocity_pseudotime'],adata_total.obs['velocity_pseudotime']) # 0.27762895
np.corrcoef(adata_total.obs['velocity_pseudotime'],adata_split2_seed320.obs['velocity_pseudotime']) # 0.81290299


out_fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/scvelo/Jun13_pseudotime/"
#adata_raw = sc.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad")
#adata_raw.uns['celltype_colors']
cell_types = adata_total.obs['celltype']
colors = dict(zip(['Blood progenitors 1', 'Blood progenitors 2', 'Erythroid1', 'Erythroid2', 'Erythroid3'],
                  ['#f9decf', '#c9a997', '#C72228', '#f79083', '#EF4E22']))

def ptime_scatter_plot(s1,s2,seed,method,name,xlab,ylab,data="ery"):
    df = pd.DataFrame({'split1':s1.obs['velocity_pseudotime'], 'split2':s2.obs['velocity_pseudotime'],
                       'cell_types':cell_types})
    corr = np.round(np.corrcoef(s1.obs['velocity_pseudotime'],s2.obs['velocity_pseudotime']),3)[0,1]
    print(corr)
    for category, color in colors.items(): plt.scatter([], [], color=color, label=category)
    plt.scatter(df['split1'], df['split2'], c=df['cell_types'].map(colors))
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title('Pseudotime '+name+' (seed='+str(seed)+', corr='+str(corr)+')')
    plt.savefig(out_fig_folder+data+"_"+method+"_seed"+str(seed)+"_pseudotime"+name+".png")
    plt.close()

ptime_scatter_plot(s1=adata_split1_seed317,s2=adata_split2_seed317,seed=317,
                   method="scv", name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=adata_split1_seed317,s2=adata_total,seed=317,
                   method="scv", name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=adata_split2_seed317,s2=adata_total,seed=317,
                   method="scv", name="split2vstotal",xlab="split2",ylab="total")


ptime_scatter_plot(s1=adata_split1_seed320,s2=adata_split2_seed320,seed=320,
                   method="scv", name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=adata_split1_seed320,s2=adata_total,seed=320,
                   method="scv", name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=adata_split2_seed320,s2=adata_total,seed=320,
                   method="scv", name="split2vstotal",xlab="split2",ylab="total")

############################################################
############################################################
## UniTVelo
out_fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/Jun13_pseudotime/"
data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/"

total_res = scv.read(data_folder+"ery_utvgenes_total.h5ad")
s1_res317 = scv.read(data_folder+"ery_utvgenes_seed317_split1.h5ad")
s2_res317 = scv.read(data_folder+"ery_utvgenes_seed317_split2.h5ad")
s1_res320 = scv.read(data_folder+"ery_utvgenes_seed320_split1.h5ad")
s2_res320 = scv.read(data_folder+"ery_utvgenes_seed320_split2.h5ad")

ptime_scatter_plot(s1=s1_res317,s2=s2_res317,seed=317,method="utv", name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=s1_res317,s2=total_res,seed=317,method="utv", name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=s2_res317,s2=total_res,seed=317,method="utv", name="split2vstotal",xlab="split2",ylab="total")

ptime_scatter_plot(s1=s1_res320,s2=s2_res320,seed=320,method="utv", name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=s1_res320,s2=total_res,seed=320,method="utv", name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=s2_res320,s2=total_res,seed=320,method="utv", name="split2vstotal",xlab="split2",ylab="total")


np.corrcoef(s1_res317.obs['velocity_pseudotime'], s2_res317.obs['velocity_pseudotime'])[0,1]
### 0.7681336936467976
np.corrcoef(s1_res320.obs['velocity_pseudotime'], s2_res320.obs['velocity_pseudotime'])[0,1]
### 0.8153660712975266

np.corrcoef(total_res.obs['velocity_pseudotime'], s1_res317.obs['velocity_pseudotime'])[0,1] ### 0.7370866473721769
np.corrcoef(total_res.obs['velocity_pseudotime'], s2_res317.obs['velocity_pseudotime'])[0,1] ### 0.9191245138388167
np.corrcoef(total_res.obs['velocity_pseudotime'], s1_res320.obs['velocity_pseudotime'])[0,1] ### 0.9196157953316147
np.corrcoef(total_res.obs['velocity_pseudotime'], s2_res320.obs['velocity_pseudotime'])[0,1] ### 0.9178693353597748







