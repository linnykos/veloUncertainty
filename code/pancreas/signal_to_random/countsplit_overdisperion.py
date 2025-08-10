import numpy as np
from scipy import sparse
from scipy.io import mmwrite
import scanpy as sc

data_folder = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/"
adata = sc.read_h5ad(data_folder+"Pancreas/endocrinogenesis_day15.h5ad")

mmwrite(data_folder+"v4_pancreas/pan_spliced.mtx", adata.layers['spliced'])
mmwrite(data_folder+"v4_pancreas/pan_unspliced.mtx", adata.layers['unspliced'])