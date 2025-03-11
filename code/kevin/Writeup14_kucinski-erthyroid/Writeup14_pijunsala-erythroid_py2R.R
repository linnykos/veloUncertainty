rm(list=ls())

library(Seurat)
library(SeuratDisk)

devtools::session_info()

data_folder <- "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup14/"
  
SeuratDisk::Convert(paste0(data_folder, "pijunsala-erthyroid_cleaned.h5ad"), 
                    dest = "h5seurat", 
                    overwrite = TRUE)

seurat_obj <- SeuratDisk::LoadH5Seurat(paste0(data_folder, "pijunsala-erthyroid_cleaned.h5seurat"))

######

# put the raw counts back in
# Load the Matrix Market format
count_mat <- Matrix::readMM(paste0(data_folder, "pijunsala-erthyroid_adata_layer_spliced.mtx"))
# Convert dgTMatrix to dgCMatrix
count_mat <- as(count_mat, "CsparseMatrix")
rownames(count_mat) <- Seurat::Cells(seurat_obj)
colnames(count_mat) <- SeuratObject::Features(seurat_obj)
count_mat <- Matrix::t(count_mat)
seurat_obj[["spliced"]] <- Seurat::CreateAssayObject(counts = count_mat)

# put the raw counts back in
# Load the Matrix Market format
count_mat <- Matrix::readMM(paste0(data_folder, "pijunsala-erthyroid_adata_layer_unspliced.mtx"))
# Convert dgTMatrix to dgCMatrix
count_mat <- as(count_mat, "CsparseMatrix")
rownames(count_mat) <- Seurat::Cells(seurat_obj)
colnames(count_mat) <- SeuratObject::Features(seurat_obj)
count_mat <- Matrix::t(count_mat)
seurat_obj[["unspliced"]] <- Seurat::CreateAssayObject(counts = count_mat)

######

# put the cell metadata back in
metadata <- read.csv(paste0(data_folder, "pijunsala-erthyroid_adata_obs_backup.csv"))
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

# convert characters into factors
for(i in 1:ncol(metadata)){
  if(is.character(metadata[,i]) && length(unique(metadata[,i])) < nrow(metadata)/2){
    metadata[,i] <- factor(metadata[,i])
  }
}
seurat_obj@meta.data <- metadata

# put the feature metadata back in
var.metadata <- read.csv(paste0(data_folder, "pijunsala-erthyroid_adata_var_backup.csv"))
rownames(var.metadata) <- var.metadata[,1]
var.metadata <- var.metadata[,-1]
seurat_obj[["RNA"]]@misc <- var.metadata

######

# put misc back in
library(jsonlite)

# Load JSON metadata
uns_metadata <- fromJSON(paste0(data_folder, "pijunsala-erthyroid_adata_metadata.json"))

# Convert character vectors back to named lists (for colors, etc.)
for (key in names(uns_metadata)) {
  if (is.character(uns_metadata[[key]])) {
    uns_metadata[[key]] <- unname(as.list(uns_metadata[[key]]))
  }
}

# Assume 'seurat_obj' is your existing Seurat object
seurat_obj@misc <- uns_metadata  # Attach metadata to @misc slot

# convert the "colors" into vectors
color_idx <- grep("_colors", names(seurat_obj@misc))
if(length(color_idx) > 0){
  for(i in color_idx){
    colname_str <- strsplit(names(seurat_obj@misc)[i], split = "_colors")[[1]][1]
    colname_str <- gsub(pattern = " ", replacement = ".", x = colname_str)
    colname_str <- gsub(pattern = "-", replacement = ".", x = colname_str)
    vec <- unlist(seurat_obj@misc[[i]])
    names(vec) <- levels(seurat_obj@meta.data[,colname_str])
    seurat_obj@misc[[i]] <- vec
  }
}

#######

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(seurat_obj,
     date_of_run, session_info,
     file = paste0(data_folder, "pijunsala-erthyroid.RData"))
