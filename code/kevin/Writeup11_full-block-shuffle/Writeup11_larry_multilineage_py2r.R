rm(list=ls()); gc(TRUE)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(biomaRt)
options(Seurat.object.assay.version = "v5")

# from https://github.com/theislab/zellkonverter/issues/38
adata <- zellkonverter::readH5AD(
  file = "~/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry_multilineage.h5ad"
)

names(adata@assays)

tmp <- Seurat::as.Seurat(adata, counts = "X", data = "X")
# converted most of the things. It creates an assay by default called "originalexp"
# see also https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html

seurat_obj <- Seurat::CreateSeuratObject(
  counts = tmp[["originalexp"]],
  data = tmp[["originalexp"]],
  meta.data = tmp@meta.data
)

seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

# put in the assays
name_vec <- names(adata@assays)
gene_vec <- SeuratObject::Features(seurat_obj[["RNA"]])

# put in the gene metafeatures
gene_metadata <- SingleCellExperiment::rowData(adata)
seurat_obj[["RNA"]]@misc <- as.data.frame(gene_metadata)

# put in the dimension reductions
name_vec <- SingleCellExperiment::reducedDimNames(adata)
for(name_val in name_vec){
  mat <- SingleCellExperiment::reducedDim(adata, name_val)
  name_val2 <- paste0("python_", name_val)
  colnames(mat) <- paste0(name_val2, "_", 1:ncol(mat))
  
  seurat_obj[[name_val2]] <- Seurat::CreateDimReducObject(embeddings = mat,
                                                          assay = "RNA")
}

# put in the metadata
metadata_list <- adata@metadata
graph_list <- numeric(0)
if(length(metadata_list) > 0){
  idx <- which(sapply(1:length(metadata_list), function(i){
    any(class(metadata_list[[i]]) %in% c("dgCMatrix", "dgRMatrix"))
  }))
  if(length(idx) > 0){
    graph_list <- metadata_list[idx]
    metadata_list <- metadata_list[-idx]
  } else {
    graph_list <- numeric(0)
  }
  seurat_obj@misc <- metadata_list
}

if(length(graph_list) > 0){
  for(name_val in names(graph_list)){
    print(paste0("Putting in graph ", name_val))
    
    seurat_obj@graphs[[name_val]] <- graph_list[[name_val]]
    rownames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)
    colnames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)
    
  }
}

Seurat::DefaultAssay(seurat_obj) <- "RNA"

# clearing environment
ls_vec <- ls()
ls_vec <- ls_vec[ls_vec != "seurat_obj"]
rm(list = ls_vec)
gc(TRUE)

################

set.seed(10)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                           selection.method = "vst",
                                           nfeatures = 2000)
seurat_obj <- Seurat::ScaleData(seurat_obj)

set.seed(10)
seurat_obj <- Seurat::RunPCA(seurat_obj,
                             assay = "RNA",
                             features = Seurat::VariableFeatures(seurat_obj),
                             verbose = FALSE)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              dims = 1:30)

####################

pdf(paste0("~/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/fig/kevin/Writeup11/Writeup11_larry_multilineage_umap.pdf"),
    onefile = T, width = 7, height = 5)

scCustomize::DimPlot_scCustom(seurat_obj,
                              group.by = "state_info",
                              reduction = "python_X_emb")

scCustomize::DimPlot_scCustom(seurat_obj,
                              group.by = "time_info",
                              reduction = "python_X_emb")

scCustomize::DimPlot_scCustom(seurat_obj,
                              group.by = "state_info",
                              reduction = "umap")

scCustomize::DimPlot_scCustom(seurat_obj,
                              group.by = "time_info",
                              reduction = "umap")

graphics.off()

################

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(seurat_obj,
     date_of_run, session_info,
     file = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/Writeup11_larry_multilineage.RData")

print("Done! :)")
