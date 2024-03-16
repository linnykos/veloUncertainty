rm(list=ls())
library(Seurat)
library(SeuratDisk)
load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/zhaoheng/out/genePathwayManifold/velocity_conversions/scvelo_pancreas_seurat.RData")

spliced_mat <- SeuratObject::LayerData(scvelo_seurat, layer = "counts", assay = "spliced")
unspliced_mat <- SeuratObject::LayerData(scvelo_seurat, layer = "counts", assay = "unspliced")

dim(spliced_mat)
dim(unspliced_mat)

spliced_mat[1:5,1:5]
unspliced_mat[1:5,1:5]

# sparsity
length(spliced_mat@x)/prod(dim(spliced_mat))
length(unspliced_mat@x)/prod(dim(unspliced_mat))

# compute the total counts
total_mat <- spliced_mat + unspliced_mat

# compute the overdispersion
vec <- as.numeric(total_mat[1,])
hist(vec, breaks = seq(min(vec), max(vec), by = 1))
df <- data.frame(x = vec)
nb_res <- MASS::glm.nb(x ~ 1, data = df)
nb_res$theta

# compute the mean and variances
mean(vec); var(vec)
mean(vec) + mean(vec)^2/nb_res$theta

# convert to AnnData via https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
seurat_simple <- Seurat::CreateSeuratObject(counts = total_mat)
seurat_simple[["RNA"]] <- as(object = seurat_simple[["RNA"]], Class = "Assay")
seurat_simple[["spliced"]] <- Seurat::CreateAssayObject(counts = spliced_mat, Class = "Assay")
seurat_simple[["unspliced"]] <- Seurat::CreateAssayObject(counts = unspliced_mat, Class = "Assay")
seurat_simple@meta.data <- scvelo_seurat@meta.data[,c("clusters_coarse", "clusters")]

SeuratDisk::SaveH5Seurat(seurat_simple,
                         filename = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/zhaoheng/out/genePathwayManifold/velocity_conversions/scvelo_pancreas_seurat_example.h5Seurat",
                         overwrite = TRUE)
SeuratDisk::Convert("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/zhaoheng/out/genePathwayManifold/velocity_conversions/scvelo_pancreas_seurat_example.h5Seurat",
                    dest = "h5ad",
                    overwrite = TRUE)

# https://www.biostars.org/p/9580132/#9580145
# https://github.com/cellgeni/sceasy
# sceasy::convertFormat(seurat_simple,
#                       from = "seurat",
#                       to = "anndata",
#                       outFile = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/zhaoheng/out/genePathwayManifold/velocity_conversions/scvelo_pancreas_seurat_example.h5ad")

