library(Seurat)
library(SeuratDisk)
library(countsplit)
library(zellkonverter)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v5")

scvelo_sc <- zellkonverter::readH5AD(file="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/ery_utv_preprocess.h5ad")

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

# now do the usual seurat processing
scvelo_seurat <- Seurat::ScaleData(scvelo_seurat)
scvelo_seurat <- Seurat::RunPCA(scvelo_seurat, features = Seurat::VariableFeatures(object = scvelo_seurat), verbose = F)
scvelo_seurat <- Seurat::RunUMAP(scvelo_seurat, dims = 1:30)

#save(scvelo_seurat,file="/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/erythroid_split/adata_erythroid.RData")

print("Object ready!")

########################
## get spliced and unspliced matrices, and a matrix of total counts
spliced_mat <- SeuratObject::LayerData(scvelo_seurat, assay = "spliced_original") 
unspliced_mat <- SeuratObject::LayerData(scvelo_seurat, assay = "unspliced_original")

print(paste0("dimension of spliced matrix is ",dim(spliced_mat)[1],"x",dim(spliced_mat)[2]," - transpose then!"))
spliced_mat <- Matrix::t(spliced_mat)
unspliced_mat <- Matrix::t(unspliced_mat)

# compute the total counts
total_mat <- spliced_mat + unspliced_mat

## Compute overdispersion parameter ####
## glm.nb(u~1) to compute the p overdispersion parameters (1 for each gene)
nb_res <- apply(total_mat, MARGIN=2, function(u) { MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(u)))$theta })

## Run countsplit on spliced and unspliced matriced ####
### helper functions
seed_split <- NULL
# spliced_split <- countsplit(spliced_mat, folds=2, epsilon=c(.5,.5), overdisps=nb_res)
get_splits <- function(mat, seed_arg, overdisps=nb_res) {
  set.seed(seed_arg)
  countsplit(mat, folds=2, epsilon=c(.5,.5), overdisps=overdisps)
}
### on Bayes
convert_object_erythroid <- function(s_mat,u_mat,fname,seed_arg) {
  print(paste0("Converting object, will transpose count matrices again!"))
  
  s_mat <- Matrix::t(s_mat)
  u_mat <- Matrix::t(u_mat)
  mat <- s_mat+u_mat
  print(paste0("dimension of the matrix stored in h5seurat and h5ad: ",dim(mat)[1],"x",dim(mat)[2]))
  
  mat <- s_mat+u_mat
  seurat_res <- Seurat::CreateSeuratObject(counts = mat)
  seurat_res[["RNA"]] <- as(object = seurat_res[["RNA"]], Class = "Assay")
  seurat_res[["spliced"]] <- Seurat::CreateAssayObject(counts = s_mat, Class = "Assay")
  seurat_res[["unspliced"]] <- Seurat::CreateAssayObject(counts = u_mat, Class = "Assay")
  seurat_res@meta.data <- data.frame("celltype"=as.character(scvelo_seurat@meta.data$celltype),
                                     "sequencing.batch"=as.character(scvelo_seurat@meta.data$sequencing.batch), ## original 1,2,3
                                     "sample"=as.character(scvelo_seurat@meta.data$sample)) ## original 2,...,37
  path <- paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/","erythroid_seed",seed_arg,fname,".h5Seurat")
  print(path)
  SeuratDisk::SaveH5Seurat(seurat_res, filename=path, overwrite = TRUE)
  SeuratDisk::Convert(path, dest = "h5ad", overwrite = TRUE)
}

### seed=317 ####
seed_split <- 317
print("Set seed = 317!")

## countsplit
spliced_split <- get_splits(spliced_mat,seed_split)
spliced_split1 <- spliced_split[[1]]
spliced_split2 <- spliced_split[[2]]

unspliced_split <- get_splits(unspliced_mat,seed_split)
unspliced_split1 <- unspliced_split[[1]]
unspliced_split2 <- unspliced_split[[2]]

## Convert object 
convert_object_erythroid(spliced_split1,unspliced_split1, "_split1_seurat",seed_split)
convert_object_erythroid(spliced_split2, unspliced_split2, "_split2_seurat",seed_split)
convert_object_erythroid(spliced_mat, unspliced_mat, "_total_seurat",seed_split)

### seed=320 ####
seed_split <- 320
print("Set seed = 320!")

## countsplit
spliced_split <- get_splits(spliced_mat,seed_split)
spliced_split1 <- spliced_split[[1]]
spliced_split2 <- spliced_split[[2]]

unspliced_split <- get_splits(unspliced_mat,seed_split)
unspliced_split1 <- unspliced_split[[1]]
unspliced_split2 <- unspliced_split[[2]]

## Convert object 
convert_object_erythroid(spliced_split1,unspliced_split1, "_split1_seurat",seed_split)
convert_object_erythroid(spliced_split2, unspliced_split2, "_split2_seurat",seed_split)
convert_object_erythroid(spliced_mat, unspliced_mat, "_total_seurat",seed_split)


##########################################
### correlation plot
s1_seed317 <- LoadH5Seurat("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/erythroid_seed317_split1_seurat.h5Seurat")
s2_seed317 <- LoadH5Seurat("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/erythroid_seed317_split2_seurat.h5Seurat")
s1_seed320 <- LoadH5Seurat("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/erythroid_seed320_split1_seurat.h5Seurat")
s2_seed320 <- LoadH5Seurat("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/ery_utv_utvgenes/erythroid_seed320_split2_seurat.h5Seurat")

spliced_seed317 <- list(SeuratObject::LayerData(s1_seed317,assay="spliced"),
                        SeuratObject::LayerData(s2_seed317,assay="spliced"))
unspliced_seed317 <- list(SeuratObject::LayerData(s1_seed317,assay="unspliced"),
                          SeuratObject::LayerData(s2_seed317,assay="unspliced"))
spliced_seed320 <- list(SeuratObject::LayerData(s1_seed320,assay="spliced"),
                        SeuratObject::LayerData(s2_seed320,assay="spliced"))
unspliced_seed320 <- list(SeuratObject::LayerData(s1_seed320,assay="unspliced"),
                          SeuratObject::LayerData(s2_seed320,assay="unspliced"))

### produce scatterplots of correlations
total = spliced_seed317[[1]] + spliced_seed317[[2]]
p <- nrow(total)
cor_spliced_seed317 <- sapply(1:nrow(spliced_seed317[[1]]),function(i){ 
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(spliced_seed317[[1]][i,]+1), log10(spliced_seed317[[2]][i,]+1) ) } )

cor_unspliced_seed317 <- sapply(1:nrow(unspliced_seed317[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(unspliced_seed317[[1]][i,]+1), log10(unspliced_seed317[[2]][i,]+1) ) } )

cor_spliced_seed320 <- sapply(1:nrow(spliced_seed320[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(spliced_seed320[[1]][i,]+1), log10(spliced_seed320[[2]][i,]+1) ) } )

cor_unspliced_seed320 <- sapply(1:nrow(unspliced_seed320[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(unspliced_seed320[[1]][i,]+1), log10(unspliced_seed320[[2]][i,]+1) ) } )

nb_res <- apply(total, MARGIN=2, function(u) { MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(u)))$theta })

zero_fraction = sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  length(which(total[j,] == 0))/ncol(total)
})
break_values <- seq(0, 1, length.out = 20)
color_palette <- grDevices::colorRampPalette(c("lightgray", "coral"))(20)
color_values <- sapply(zero_fraction, function(val){
  color_palette[which.min(abs(val - break_values))]
})

png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/corr_utvgenes_seed317.png"), height = 1000, width = 1800, units = "px", res = 300)
par(mfrow=c(1,2))
plot(cor_spliced_seed317, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Spliced",ylab="corr",cex=.7,cex.main=.8) 
plot(cor_unspliced_seed317, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Unspliced",ylab="corr",cex=.7,cex.main=.8)
graphics.off()

png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/erythroid/unitvelo_utvgenes/corr_utvgenes_seed320.png"), height = 1000, width = 1800, units = "px", res = 300)
par(mfrow=c(1,2))
plot(cor_spliced_seed320, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Spliced",ylab="corr",cex=.7,cex.main=.8) 
plot(cor_unspliced_seed320, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Unspliced",ylab="corr",cex=.7,cex.main=.8)
graphics.off()





