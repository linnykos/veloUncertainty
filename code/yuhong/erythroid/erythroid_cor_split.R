### run on local R
library(Seurat)
library(SeuratDisk)
library(SeuratObject)

setwd("/Users/wakeup/Downloads/UW_23-25/proj_RNA/data/erythroid")

split1_seed317 <- LoadH5Seurat("./split_data/scvelo_erythroid_split1_seurat_seed317.h5Seurat")
split2_seed317 <- LoadH5Seurat("./split_data/scvelo_erythroid_split2_seurat_seed317.h5Seurat")
split1_seed320 <- LoadH5Seurat("./split_data/scvelo_erythroid_split1_seurat_seed320.h5Seurat")
split2_seed320 <- LoadH5Seurat("./split_data/scvelo_erythroid_split2_seurat_seed320.h5Seurat")

spliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="spliced"),
                        SeuratObject::LayerData(split2_seed317,assay="spliced"))
unspliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed317,assay="unspliced"))
spliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="spliced"),
                        SeuratObject::LayerData(split2_seed320,assay="spliced"))
unspliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed320,assay="unspliced"))

### produce scatterplots of correlations
cor_spliced_seed317 <- sapply(1:nrow(spliced_seed317[[1]]),function(i){ 
  cor( log10(spliced_seed317[[1]][i,]+1), log10(spliced_seed317[[2]][i,]+1) ) } )

cor_unspliced_seed317 <- sapply(1:nrow(unspliced_seed317[[1]]),function(i){
  cor( log10(unspliced_seed317[[1]][i,]+1), log10(unspliced_seed317[[2]][i,]+1) ) } )

cor_spliced_seed320 <- sapply(1:nrow(spliced_seed320[[1]]),function(i){
  cor( log10(spliced_seed320[[1]][i,]+1), log10(spliced_seed320[[2]][i,]+1) ) } )

cor_unspliced_seed320 <- sapply(1:nrow(unspliced_seed320[[1]]),function(i){
  cor( log10(unspliced_seed320[[1]][i,]+1), log10(unspliced_seed320[[2]][i,]+1) ) } )

plot(cor_spliced_seed317, ylim=c(-1,1))
dev.copy(png,filename="/Users/wakeup/Downloads/UW_23-25/proj_RNA/git/veloUncertainty/fig/yuhong/erythroid/cor_spliced_seed317.png")
dev.off()

plot(cor_unspliced_seed317, ylim=c(-1,1))
dev.copy(png,filename="/Users/wakeup/Downloads/UW_23-25/proj_RNA/git/veloUncertainty/fig/yuhong/erythroid/cor_unspliced_seed317.png")
dev.off()

plot(cor_spliced_seed320, ylim=c(-1,1))
dev.copy(png,filename="/Users/wakeup/Downloads/UW_23-25/proj_RNA/git/veloUncertainty/fig/yuhong/erythroid/cor_spliced_seed320.png")
dev.off()

plot(cor_unspliced_seed320, ylim=c(-1,1))
dev.copy(png,filename="/Users/wakeup/Downloads/UW_23-25/proj_RNA/git/veloUncertainty/fig/yuhong/erythroid/cor_unspliced_seed320.png")
dev.off()


#### NA's in correlations between unspliced splits - due to (at least) one split matrix contains all zero counts for a particular gene
which(is.na(cor_unspliced_seed317))
## 541 1377 1559 1612
which(is.na(cor_unspliced_seed320))
## 541 1377 1612
check_nonzero_seed317 <- sapply(1:nrow(unspliced_seed317[[1]]),function(i){ (sum(unspliced_seed317[[1]][i,])>0 & sum(unspliced_seed317[[2]][i,])>0) })
which(!check_nonzero_seed317) ### 541 1377 1559 1612
check_nonzero_seed320 <- sapply(1:nrow(unspliced_seed320[[1]]),function(i){ (sum(unspliced_seed320[[1]][i,])>0 & sum(unspliced_seed320[[2]][i,])>0) })
which(!check_nonzero_seed320) ### 541 1377 1612

summary(cor_spliced_seed317)
summary(cor_unspliced_seed317)
summary(cor_spliced_seed320)
summary(cor_unspliced_seed320)

