import scvelo as scv
import scanpy as sc

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes/"
figure_folder="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/velo_conf/"

label='clusters'
color_map='clusters_colors'

adata = scv.read(data_folder+'pan_utv_preprocess.h5ad')

total_res = scv.read(data_folder+"pan_utvgenes_total.h5ad")
s1_res317 = scv.read(data_folder+"pan_utvgenes_seed317_split1.h5ad")
s2_res317 = scv.read(data_folder+"pan_utvgenes_seed317_split2.h5ad")
s1_res320 = scv.read(data_folder+"pan_utvgenes_seed320_split1.h5ad")
s2_res320 = scv.read(data_folder+"pan_utvgenes_seed320_split2.h5ad")

def plot_velo_conf(data,fig_name):
    scv.pl.scatter(data, color=label, cmap=color_map, basis="pca",
                   save=figure_folder+"pan_utvgenes_pca_scatter_"+fig_name+".png")
    scv.tl.velocity_confidence(data)
    scv.pl.scatter(data, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], basis="pca",
                   save=figure_folder+"pan_utvgenes_pca_velo_confidence_"+fig_name+".png")
    data.obsm['X_umap'] = adata.obsm['X_umap'].copy()
    scv.pl.scatter(data, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], basis="umap",
                   save=figure_folder+"pan_utvgenes_umapOriginal_velo_confidence_"+fig_name+".png")
    del data.obsm['X_umap']
    sc.tl.umap(data)
    scv.pl.scatter(data, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], basis="umap",
                   save=figure_folder+"pan_utvgenes_umapRecover_velo_confidence_"+fig_name+".png")
    scv.pl.scatter(data, color=label, cmap=color_map, basis="umap",
                   save=figure_folder+"pan_utvgenes_umapRecover_scatter_"+fig_name+".png")

scv.pl.scatter(adata, color=label, cmap=color_map, basis="umap",
                   save=figure_folder+"pan_utvgenes_umapOriginal_preprocess.png")
scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], basis="pca",
               save=figure_folder+"pan_utvgenes_pca_velo_confidence_preprocess.png")
scv.pl.scatter(adata, color=label, cmap=color_map, basis="pca",
               save=figure_folder+"pan_utvgenes_pca_scatter_preprocess.png")
scv.pl.scatter(adata, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], basis="umap", 
               save=figure_folder+"pan_utvgenes_umapOriginal_velo_confidence_preprocess.png")
scv.pl.scatter(adata, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], basis="umap", vmin=0, vmax=1,
               save=figure_folder+"pan_utvgenes_umapOriginal_velo_confidence_preprocess_test.png")

plot_velo_conf(total_res, "total")
plot_velo_conf(s1_res317, "seed317_s1")
plot_velo_conf(s2_res317, "seed317_s2")
plot_velo_conf(s1_res320, "seed320_s1")
plot_velo_conf(s2_res320, "seed320_s2")
