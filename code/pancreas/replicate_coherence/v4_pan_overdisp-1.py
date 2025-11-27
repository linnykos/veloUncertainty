import scanpy as sc
import scipy.io

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"

adata = sc.read(data_folder+'Pancreas/endocrinogenesis_day15.h5ad')
spliced_matrix = adata.layers['spliced']
unspliced_matrix = adata.layers['unspliced']

# Write to .mtx
scipy.io.mmwrite(data_folder+'v4_pancreas/pan_spliced.mtx', spliced_matrix)
scipy.io.mmwrite(data_folder+'v4_pancreas/pan_unspliced.mtx', unspliced_matrix)

