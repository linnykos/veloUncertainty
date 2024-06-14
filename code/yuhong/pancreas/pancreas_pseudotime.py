import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

adata_pre = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/adata_pancreas_preprocess.h5ad")

def run_scv(adata):
    scv.pp.normalize_per_cell(adata)
    scv.pp.log1p(adata)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    return adata

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/"
adata_total = scv.read(data_folder+'pancreas_seed317_total_seurat.h5ad')
run_scv(adata_total)
adata_total.write_h5ad(data_folder+'pancreas_seed317_total.h5ad')

# read split counts data
adata_split1_seed317 = scv.read(data_folder+'pancreas_seed317_split1_seurat.h5ad')
adata_split2_seed317 = scv.read(data_folder+'pancreas_seed317_split2_seurat.h5ad')
run_scv(adata_split1_seed317)
adata_split1_seed317.write_h5ad(data_folder+'pancreas_seed317_split1.h5ad')
run_scv(adata_split2_seed317)
adata_split2_seed317.write_h5ad(data_folder+'pancreas_seed317_split2.h5ad')

np.corrcoef(adata_split1_seed317.obs['velocity_pseudotime'],adata_split2_seed317.obs['velocity_pseudotime'])
### 0.28357643
np.corrcoef(adata_total.obs['velocity_pseudotime'],adata_split1_seed317.obs['velocity_pseudotime']) # 0.2462441
np.corrcoef(adata_total.obs['velocity_pseudotime'],adata_split2_seed317.obs['velocity_pseudotime']) # 0.99438046

out_fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun13_pseudotime/"

### seed=320
adata_split1_seed320 = scv.read(data_folder+'pancreas_seed320_split1_seurat.h5ad')
adata_split2_seed320 = scv.read(data_folder+'pancreas_seed320_split2_seurat.h5ad')
run_scv(adata_split1_seed320)
run_scv(adata_split2_seed320)
adata_split1_seed320.write_h5ad(data_folder+'pancreas_seed320_split1.h5ad')
adata_split2_seed320.write_h5ad(data_folder+'pancreas_seed320_split2.h5ad')

np.corrcoef(adata_split1_seed320.obs['velocity_pseudotime'],adata_split2_seed320.obs['velocity_pseudotime'])
### 0.70963656
np.corrcoef(adata_split1_seed320.obs['velocity_pseudotime'],adata_total.obs['velocity_pseudotime']) # 0.99084514
np.corrcoef(adata_total.obs['velocity_pseudotime'],adata_split2_seed320.obs['velocity_pseudotime']) # 0.66869244

cell_types = adata_total.obs['clusters']
colors = dict(zip(['Ductal', 'Ngn3 low EP', 'Ngn3 high EP', 'Pre-endocrine', 'Beta','Alpha', 'Delta', 'Epsilon'],
                  ['#8fbc8f', '#f4a460', '#fdbf6f', '#ff7f00', '#b2df8a', '#1f78b4','#6a3d9a', '#cab2d6']))

def ptime_scatter_plot(s1,s2,seed,method,name,xlab,ylab,data="pan"):
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
# seed=317
ptime_scatter_plot(s1=adata_split1_seed317, s2=adata_split2_seed317, seed=317, 
                   method="scv", name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=adata_split1_seed317, s2=adata_total, seed=317, 
                   method="scv", name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=adata_split2_seed317, s2=adata_total, seed=317, 
                   method="scv", name="split2vstotal",xlab="split2",ylab="total")

# seed=320
ptime_scatter_plot(s1=adata_split1_seed320, s2=adata_split2_seed320, seed=320, 
                   method="scv", name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=adata_split1_seed320, s2=adata_total, seed=320, 
                   method="scv", name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=adata_split2_seed320, s2=adata_total, seed=320, 
                   method="scv", name="split2vstotal",xlab="split2",ylab="total")


np.corrcoef(adata_split1_seed317.obs['velocity_self_transition'],adata_split2_seed317.obs['velocity_self_transition'])[0,1]
### 0.38179338658579187
np.corrcoef(adata_split1_seed320.obs['velocity_self_transition'],adata_split2_seed320.obs['velocity_self_transition'])[0,1]
### 0.3570479000695066

np.corrcoef(adata_total.obs['velocity_self_transition'],adata_split1_seed317.obs['velocity_self_transition'])[0,1] # 0.48407347381141697
np.corrcoef(adata_total.obs['velocity_self_transition'],adata_split2_seed317.obs['velocity_self_transition'])[0,1] # 0.4844992911983179

np.corrcoef(adata_total.obs['velocity_self_transition'],adata_split1_seed320.obs['velocity_self_transition'])[0,1] # 0.5599229116344328
np.corrcoef(adata_total.obs['velocity_self_transition'],adata_split2_seed320.obs['velocity_self_transition'])[0,1] # 0.4737363317552139


############################################################
############################################################
## UniTVelo
out_fig_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/Jun13_pseudotime/"

adata_total_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_total.h5ad') # 3696 × 1945
split1_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_split1.h5ad') # 3696 × 734
split2_seed317_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed317_split2.h5ad') # 3696 × 731
split1_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split1.h5ad') 
split2_seed320_res = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/pan_utvgenes_seed320_split2.h5ad')

## cell-specific pseudotime: split1_seed317_res.obs['velocity_pseudotime']
np.corrcoef(split1_seed317_res.obs['velocity_pseudotime'], split2_seed317_res.obs['velocity_pseudotime'])[0,1]
### 0.9942812176441053
np.corrcoef(split1_seed320_res.obs['velocity_pseudotime'], split2_seed320_res.obs['velocity_pseudotime'])[0,1]
### 0.9949345690567206
ptime_scatter_plot(s1=split1_seed317_res,s2=split2_seed317_res,seed=317,method="utv",
                   name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=split1_seed317_res,s2=adata_total_res,seed=317,method="utv",
                   name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=split2_seed317_res,s2=adata_total_res,seed=317,method="utv",
                   name="split2vstotal",xlab="split2",ylab="total")

ptime_scatter_plot(s1=split1_seed320_res,s2=split2_seed320_res,seed=320,method="utv",
                   name="split1vs2",xlab="split1",ylab="split2")
ptime_scatter_plot(s1=split1_seed320_res,s2=adata_total_res,seed=320,method="utv",
                   name="split1vstotal",xlab="split1",ylab="total")
ptime_scatter_plot(s1=split2_seed320_res,s2=adata_total_res,seed=320,method="utv",
                   name="split2vstotal",xlab="split2",ylab="total")

## split vs total
np.corrcoef(split1_seed317_res.obs['velocity_pseudotime'], adata_total_res.obs['velocity_pseudotime'])[0,1] # 0.9846630025417334
np.corrcoef(adata_total_res.obs['velocity_pseudotime'], split2_seed317_res.obs['velocity_pseudotime'])[0,1] # 0.9844597712051569

np.corrcoef(split1_seed320_res.obs['velocity_pseudotime'], adata_total_res.obs['velocity_pseudotime'])[0,1] # 0.9817933073700154
np.corrcoef(adata_total_res.obs['velocity_pseudotime'], split2_seed320_res.obs['velocity_pseudotime'])[0,1] # 0.9782390150780548



