rm(list=ls())
library(Seurat)
load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/zhaoheng/out/genePathwayManifold/velocity_conversions/scvelo_pancreas_seurat.RData")

spliced_mat <- SeuratObject::LayerData(scvelo_seurat, layer = "data", assay = "spliced_original")
unspliced_mat <- SeuratObject::LayerData(scvelo_seurat, layer = "data", assay = "spliced_original")
