import scvelo as scv
import unitvelo as utv
import scanpy as sc
import tf_keras
import os
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt


### did not run pp.neighbors
# the below script uses the environment: "utvClone"

velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = False
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.ASSIGN_POS_U = True
os.environ["TF_USE_LEGACY_KERAS"]="1"
# https://github.com/tensorflow/tensorflow/issues/62337
# To solve the following error: 
## ImportError: `keras.optimizers.legacy` is not supported in Keras 3. 
## When using `tf.keras`, to continue using a `tf.keras.optimizers.legacy` optimizer, 
## you can install the `tf_keras` package (Keras 2) and set the environment variable `TF_USE_LEGACY_KERAS=True` 
## to configure TensorFlow to use `tf_keras` when accessing `tf.keras`.
# Tried this but did not work:
## velo_config.TF_USE_LEGACY_KERAS=True

label='clusters'
adata = scv.read("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/simulation/sim3/sim3good.h5ad")
scv.pp.filter_genes(adata, min_shared_counts=20) 
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) 

adata.var['highly_variable'] = True
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# To solve the error: ValueError: '_index' is a reserved name for dataframe columns.
# ref: https://github.com/theislab/scvelo/issues/255
path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup8_simulation/tmp/tmp_sim3good.h5ad"
adata.write_h5ad(filename=path)
adata_res = utv.run_model(path, label, config_file=velo_config)
#sc.tl.umap(adata_res)
scv.pl.velocity_embedding_stream(adata_res,basis="umap",color="true_t",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/simulation/sim3/sim3good_utv.png")

split1_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed317_split1_seurat.h5ad')
split1_seed317.var['highly_variable'] = adata.var['highly_variable'].copy()
split1_seed317.__dict__['_raw'].__dict__['_var'] = split1_seed317.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split1_seed317_path = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup2_utv/tmp/scvelo_pancreas_split1_seurat_seed317.h5ad'
split1_seed317.write_h5ad(filename=split1_seed317_path)
split1_seed317_res = utv.run_model(split1_seed317_path, label, config_file=velo_config)
sc.tl.umap(split1_seed317_res)
scv.pl.velocity_embedding_stream(split1_seed317_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_split1.png")
print("seed317 split1 done!")

split2_seed317 = scv.read('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split/pancreas_seed317_split2_seurat.h5ad')
split2_seed317.var['highly_variable'] = adata.var['highly_variable'].copy()
split2_seed317.__dict__['_raw'].__dict__['_var'] = split2_seed317.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
split2_seed317_path = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup2_utv/tmp/scvelo_pancreas_split2_seurat_seed317.h5ad'
split2_seed317.write_h5ad(filename=split2_seed317_path)
split2_seed317_res = utv.run_model(split2_seed317_path, label, config_file=velo_config)
sc.tl.umap(split2_seed317_res)
scv.pl.velocity_embedding_stream(split2_seed317_res,basis="umap",color="clusters",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_split2.png")
print("seed317 split2 done!")

#split1_seed317_res.layers["velocity_rmNA"] = np.nan_to_num(split1_seed317_res.layers['velocity'], nan=0)
#split2_seed317_res.layers["velocity_rmNA"] = np.nan_to_num(split2_seed317_res.layers['velocity'], nan=0)
common_genes_317 = np.intersect1d(split1_seed317_res.var['features'], split2_seed317_res.var['features'])
#indices_seed317 = split1_seed317_res.var['features'][split1_seed317_res.var['features'].isin(set(common_genes_317))].index.tolist()
Ngenes_317s1 = len(split1_seed317_res.var['features'])
Ngenes_317s2 = len(split2_seed317_res.var['features'])
Ngenes_317common = len(common_genes_317)

# Here assuming order of genes in layers['velocity'] is the same as in var['features']
### split1_seed317_res.layers['velocity'].shape = (3696, 734)
### split2_seed317_res.layers['velocity'].shape = (3696, 731)
df1 = pd.DataFrame(split1_seed317_res.layers['velocity'], columns=split1_seed317_res.var['features'].tolist())
df2 = pd.DataFrame(split2_seed317_res.layers['velocity'], columns=split2_seed317_res.var['features'].tolist())
cos_sim_seed317 = np.diag(cosine_similarity(df1[common_genes_317],df2[common_genes_317]))
#cos_sim_seed317 = np.diag(cosine_similarity(df1[indices_seed317],df2[indices_seed317]))

# Create histogram
plt.clf()
plt.hist(cos_sim_seed317, bins=30, edgecolor='black')  # Adjust bins and edgecolor as needed
## add mean
mean_seed317 = np.mean(cos_sim_seed317)
plt.axvline(mean_seed317, color='red', linestyle='dashed', linewidth=1)
## add number of genes used in each split
plt.text(.2, 400, 'mean cosine similarity = '+str(mean_seed317), color='blue', fontsize=10)
plt.text(.2, 350, 'split1 number of genes = '+str(Ngenes_317s1), color='blue', fontsize=10)
plt.text(.2, 300, 'split2 number of genes = '+str(Ngenes_317s2), color='blue', fontsize=10)
## add labels and title
plt.xlabel('cosine similarity (seed317)')
plt.ylabel('Frequency')
plt.title('Histogram of cosine similarity, pan+utv_scvgenes, Ngenes='+str(Ngenes_317common))
plt.savefig('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_cos_similarity_hist.png')
plt.clf()

## Plot on UMAP
adata_total_res.obs['cos_sim_317'] = cos_sim_seed317
adata_total_res.obs['cos_sim_317'] = pd.DataFrame(adata_total_res.obs['cos_sim_317'])
scv.pl.velocity_embedding_stream(adata_total_res,basis="umap",color="cos_sim_317",cmap="coolwarm",
                                 save="/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo/unitvelo_seed317_cos_similarity.png")
del adata_total_res.obs['cos_sim_317']

