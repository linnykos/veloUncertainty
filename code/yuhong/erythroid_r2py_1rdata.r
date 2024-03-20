library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v5")

scvelo_sc <- zellkonverter::readH5AD(file="./adata_erythroid_v0.h5ad")

print("Dataset of erythroid loaded!")
# names(scvelo_sc@assays)
### [1] "X"                  "Ms"                 "Mu"                
### [4] "spliced"            "spliced_original"   "unspliced"         
### [7] "unspliced_original" "variance_velocity"  "velocity"          

tmp <- Seurat::as.Seurat(scvelo_sc, counts = "X", data = "X")

scvelo_seurat <- Seurat::CreateSeuratObject(counts = tmp[["originalexp"]], data = tmp[["originalexp"]], meta.data = tmp@meta.data)

scvelo_seurat[["RNA"]] <- as(object = scvelo_seurat[["RNA"]], Class = "Assay5")

# put in the assays
name_vec <- names(scvelo_sc@assays)
gene_vec <- SeuratObject::Features(scvelo_seurat[["RNA"]])
for(i in which(name_vec != "X")){
  name_val <- name_vec[i]
  print(paste0("Working on ", name_val))
  mat <- SummarizedExperiment::assay(scvelo_sc, name_val)
  if(is.matrix(mat) && length(which(mat == 0)) > prod(dim(mat))/2) {mat <- Matrix::Matrix(mat, sparse = T)}
  if(length(colnames(mat)) == 0) {colnames(mat) <- colnames(scvelo_seurat)}
  if(length(rownames(mat)) == 0) {rownames(mat) <- gene_vec}
  scvelo_seurat[[name_val]] <- SeuratObject::CreateAssay5Object(data = mat)
}
### Warning: No layers found matching search pattern provided

# replace the counts assay in RNA with the sum of spliced_original and unspliced_original
mat <- SeuratObject::LayerData(scvelo_seurat, assay = "spliced_original", layer = "data") + SeuratObject::LayerData(scvelo_seurat, assay = "unspliced_original", layer = "data")
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

  scvelo_seurat[[name_val2]] <- Seurat::CreateDimReducObject(embeddings = mat, assay = "RNA")
}
#### "velocity_umap": Warning message: Key ‘scvelo_’ taken, using ‘scveloxumap_’ instead 


# put in the metadata
metadata_list <- scvelo_sc@metadata
idx <- which(sapply(1:length(metadata_list), function(i){class(metadata_list[[i]]) %in% c("dgCMatrix", "dgRMatrix")}))
graph_list <- metadata_list[idx]
metadata_list <- metadata_list[-idx]
scvelo_seurat@misc <- metadata_list

for(name_val in names(graph_list)){
  print(paste0("Putting in graph ", name_val))
  scvelo_seurat@graphs[[name_val]] <- graph_list[[name_val]]
}

Seurat::DefaultAssay(scvelo_seurat) <- "RNA"
Seurat::VariableFeatures(scvelo_seurat) <- SeuratObject::Features(scvelo_seurat)


######### Did not run these two steps
# for the cluster color specifically, add the names for convenience
# names(scvelo_seurat@misc[["clusters_colors"]]) <- sort(unique(scvelo_seurat$clusters))
# names(scvelo_seurat@misc[["clusters_coarse_colors"]]) <- sort(unique(scvelo_seurat$clusters_coarse))

# now do the usual seurat processing
scvelo_seurat <- Seurat::ScaleData(scvelo_seurat)
scvelo_seurat <- Seurat::RunPCA(scvelo_seurat, features = Seurat::VariableFeatures(object = scvelo_seurat), verbose = F)
scvelo_seurat <- Seurat::RunUMAP(scvelo_seurat, dims = 1:30)

save(scvelo_seurat,file="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/adata_erythroid.RData")
