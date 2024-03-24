rm(list=ls())
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v5")

# from https://github.com/theislab/zellkonverter/issues/38
scvelo_sc <- zellkonverter::readH5AD(
  file = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/zhaoheng/out/genePathwayManifold/velocity_conversions/scvelo_pancreas_dataset.h5ad"
)

names(scvelo_sc@assays)

tmp <- Seurat::as.Seurat(scvelo_sc, counts = "X", data = "X")
# converted most of the things. It creates an assay by default called "originalexp"
# see also https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html

scvelo_seurat <- Seurat::CreateSeuratObject(
  counts = tmp[["originalexp"]],
  data = tmp[["originalexp"]],
  meta.data = tmp@meta.data
)

scvelo_seurat[["RNA"]] <- as(object = scvelo_seurat[["RNA"]], Class = "Assay5")

# put in the assays
name_vec <- names(scvelo_sc@assays)
gene_vec <- SeuratObject::Features(scvelo_seurat[["RNA"]])
for(i in which(name_vec != "X")){
  name_val <- name_vec[i]
  print(paste0("Working on ", name_val))

  mat <- SummarizedExperiment::assay(scvelo_sc, name_val)
  if(is.matrix(mat) && length(which(mat == 0)) > prod(dim(mat))/2) mat <- Matrix::Matrix(mat, sparse = T)
  if(length(colnames(mat)) == 0) colnames(mat) <- colnames(scvelo_seurat)
  if(length(rownames(mat)) == 0) rownames(mat) <- gene_vec
  scvelo_seurat[[name_val]] <- SeuratObject::CreateAssay5Object(data = mat)
}

# replace the counts assay in RNA with the sum of spliced_original and unspliced_original
mat <- SeuratObject::LayerData(scvelo_seurat, assay = "spliced_original", layer = "data")
SeuratObject::LayerData(scvelo_seurat, assay = "spliced", layer = "counts") <- mat
scvelo_seurat[["spliced_original"]] <- NULL

mat <- SeuratObject::LayerData(scvelo_seurat, assay = "unspliced_original", layer = "data")
SeuratObject::LayerData(scvelo_seurat, assay = "unspliced", layer = "counts") <- mat
scvelo_seurat[["unspliced_original"]] <- NULL

mat <- SeuratObject::LayerData(scvelo_seurat, assay = "spliced", layer = "counts") + SeuratObject::LayerData(scvelo_seurat, assay = "unspliced", layer = "counts")
SeuratObject::LayerData(scvelo_seurat, assay = "RNA", layer = "counts") <- mat

mat <- SeuratObject::LayerData(scvelo_seurat, assay = "spliced", layer = "data") + SeuratObject::LayerData(scvelo_seurat, assay = "unspliced", layer = "data")
SeuratObject::LayerData(scvelo_seurat, assay = "RNA", layer = "data") <- mat

# put in the gene metafeatures
gene_metadata <- SingleCellExperiment::rowData(scvelo_sc)
scvelo_seurat[["RNA"]]@misc <- as.data.frame(gene_metadata)

# put in the dimension reductions
name_vec <- SingleCellExperiment::reducedDimNames(scvelo_sc)
for(name_val in name_vec){
  mat <- SingleCellExperiment::reducedDim(scvelo_sc, name_val)
  name_val2 <- paste0("scvelo_", name_val)
  colnames(mat) <- paste0(name_val2, "_", 1:ncol(mat))

  scvelo_seurat[[name_val2]] <- Seurat::CreateDimReducObject(embeddings = mat,
                                                             assay = "RNA")
}

# put in the metadata
metadata_list <- scvelo_sc@metadata
idx <- which(sapply(1:length(metadata_list), function(i){
  class(metadata_list[[i]]) %in% c("dgCMatrix", "dgRMatrix")
}))
graph_list <- metadata_list[idx]
metadata_list <- metadata_list[-idx]
scvelo_seurat@misc <- metadata_list

for(name_val in names(graph_list)){
  print(paste0("Putting in graph ", name_val))

  scvelo_seurat@graphs[[name_val]] <- graph_list[[name_val]]
  rownames(scvelo_seurat@graphs[[name_val]]) <- SeuratObject::Cells(scvelo_seurat)
  colnames(scvelo_seurat@graphs[[name_val]]) <- SeuratObject::Cells(scvelo_seurat)

}

Seurat::DefaultAssay(scvelo_seurat) <- "RNA"
Seurat::VariableFeatures(scvelo_seurat) <- SeuratObject::Features(scvelo_seurat)

# for the cluster color specifically, add the names for convenience
names(scvelo_seurat@misc[["clusters_colors"]]) <- sort(unique(scvelo_seurat$clusters))
names(scvelo_seurat@misc[["clusters_coarse_colors"]]) <- sort(unique(scvelo_seurat$clusters_coarse))

# now do the usual seurat processing
Seurat::DefaultAssay(scvelo_seurat) <- "RNA"
scvelo_seurat <- Seurat::ScaleData(scvelo_seurat)
scvelo_seurat <- Seurat::RunPCA(scvelo_seurat,
                                features = Seurat::VariableFeatures(object = scvelo_seurat),
                                verbose = F)
scvelo_seurat <- Seurat::RunUMAP(scvelo_seurat,
                                 dims = 1:30)

####################

Seurat::DimPlot(scvelo_seurat,
                group.by = "clusters",
                cols = scvelo_seurat@misc[["clusters_colors"]],
                reduction = "scvelo_X_umap")

Seurat::DimPlot(scvelo_seurat,
                group.by = "clusters",
                cols = scvelo_seurat@misc[["clusters_colors"]],
                reduction = "umap")

Seurat::FeaturePlot(scvelo_seurat,
                    features = "velocity_pseudotime",
                    reduction = "scvelo_X_umap")

Seurat::FeaturePlot(scvelo_seurat,
                    features = "velocity_confidence",
                    reduction = "scvelo_X_umap")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Basic scVelo pancreas dataset processed following https://scvelo.readthedocs.io/en/stable/VelocityBasics.html",
              "from the script in https://github.com/zhli12/DifferentialMarkers/blob/813502e3221598cd1e361051d5bc91a087c58694/scvelo-converstion-to-R.ipynb",
              "and https://github.com/zhli12/DifferentialMarkers/blob/813502e3221598cd1e361051d5bc91a087c58694/scvelo-converstion-to-R_part2.R",
              "which uses the R packages zellkonverter and SingleCellExperiment.")

save(scvelo_seurat,
     date_of_run, session_info, note,
     file = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/zhaoheng/out/genePathwayManifold/velocity_conversions/scvelo_pancreas_seurat.RData")

print("Done! :)")
