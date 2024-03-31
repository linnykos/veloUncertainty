library(Seurat)
library(SeuratDisk)
library(SeuratObject)

setwd("~/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split")

split1_seed317 <- LoadH5Seurat("./scvelo_erythroid_split1_seurat_seed317.h5Seurat")
split2_seed317 <- LoadH5Seurat("./scvelo_erythroid_split2_seurat_seed317.h5Seurat")
split1_seed320 <- LoadH5Seurat("./scvelo_erythroid_split1_seurat_seed320.h5Seurat")
split2_seed320 <- LoadH5Seurat("./scvelo_erythroid_split2_seurat_seed320.h5Seurat")

## extract spliced and unspliced matrices, stored in list
spliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="spliced"),
                        SeuratObject::LayerData(split2_seed317,assay="spliced"))
unspliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed317,assay="unspliced"))
spliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="spliced"),
                        SeuratObject::LayerData(split2_seed320,assay="spliced"))
unspliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed320,assay="unspliced"))

save(spliced_seed317,unspliced_seed317,spliced_seed320,unspliced_seed320,
     file="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/adata_erythroid.RData")


