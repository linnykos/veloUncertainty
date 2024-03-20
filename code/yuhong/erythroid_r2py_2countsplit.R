rm(list=ls())

load("/Users/wakeup/Downloads/UW_23-25/proj_RNA/data/erythroid/adata_erythroid.RData")

library(Seurat)
library(SeuratDisk)
library(countsplit)

## get spliced and unspliced matrices, and a matrix of total counts
spliced_mat <- SeuratObject::LayerData(scvelo_seurat, assay = "spliced_original") # assay="spliced", layer="count"
unspliced_mat <- SeuratObject::LayerData(scvelo_seurat, assay = "unspliced_original")
# compute the total counts
total_mat <- spliced_mat + unspliced_mat

# dim
dim(spliced_mat) ### [1] 1949 9815
dim(unspliced_mat)
# sparsity
length(spliced_mat@x)/prod(dim(spliced_mat)) ### 0.3862658
length(unspliced_mat@x)/prod(dim(unspliced_mat)) ### 0.06445538

## Compute overdispersion parameter ####
## glm.nb(u~1) to compute the p overdispersion parameters (1 for each gene)
nb_res <- apply(total_mat, MARGIN=2, function(u) { MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(u)))$theta })

## Run countsplit on spliced and unspliced matriced ####
set.seed(317)
spliced_split <- countsplit(spliced_mat, folds=2, epsilon=c(.5,.5), overdisps=nb_res)
spliced_split1 <- spliced_split[[1]]
spliced_split2 <- spliced_split[[2]]

set.seed(317)
unspliced_split <- countsplit(unspliced_mat, folds=2, epsilon=c(.5,.5), overdisps=nb_res)
unspliced_split1 <- unspliced_split[[1]]
unspliced_split2 <- unspliced_split[[2]]

identical(spliced_split1+spliced_split2,spliced_mat)
identical(unspliced_split1+unspliced_split2,unspliced_mat)
identical(spliced_split1@p,spliced_split2@p)
identical(unspliced_split1@p,unspliced_split2@p)


###### for test purpose only 
create_split_object_erythroid <- function(s_mat,u_mat) {
  mat <- s_mat+u_mat
  seurat_res <- Seurat::CreateSeuratObject(counts = mat)
  seurat_res[["RNA"]] <- as(object = seurat_res[["RNA"]], Class = "Assay")
  seurat_res[["spliced"]] <- Seurat::CreateAssayObject(counts = s_mat, Class = "Assay")
  seurat_res[["unspliced"]] <- Seurat::CreateAssayObject(counts = u_mat, Class = "Assay")
  seurat_res@meta.data <- data.frame("celltype"=as.character(scvelo_seurat@meta.data$celltype))
  seurat_res
}
scvelo_seurat_split1 <- create_split_object_erythroid(spliced_split1,unspliced_split1)
View(scvelo_seurat_split1)
###### 

## Convert object ####
convert_object_erythroid <- function(s_mat,u_mat,fname) {
  mat <- s_mat+u_mat
  seurat_res <- Seurat::CreateSeuratObject(counts = mat)
  seurat_res[["RNA"]] <- as(object = seurat_res[["RNA"]], Class = "Assay")
  seurat_res[["spliced"]] <- Seurat::CreateAssayObject(counts = s_mat, Class = "Assay")
  seurat_res[["unspliced"]] <- Seurat::CreateAssayObject(counts = u_mat, Class = "Assay")
  seurat_res@meta.data <- data.frame("celltype"=as.character(scvelo_seurat@meta.data$celltype))
  path <- paste0("/Users/wakeup/Downloads/UW_23-25/proj_RNA/data/erythroid/",fname,".h5Seurat")
  print(path)
  SeuratDisk::SaveH5Seurat(seurat_res, filename=path, overwrite = TRUE)
  SeuratDisk::Convert(path, dest = "h5ad", overwrite = TRUE)
}

convert_object_erythroid(spliced_split1,unspliced_split1, "scvelo_erythroid_split1_seurat")
convert_object_erythroid(spliced_split2, unspliced_split2, "scvelo_erythroid_split2_seurat")
convert_object_erythroid(spliced_mat, unspliced_mat, "scvelo_erythroid_total_seurat")


## Spliced splits
### both non-zero
sum(spliced_split1@x>0 & spliced_split2@x>0)/length(spliced_split1@x)

### one of two splits non-zero and the other zero
sum(spliced_split1@x>0 & spliced_split2@x==0)/length(spliced_split1@x)
sum(spliced_split1@x==0 & spliced_split2@x>0)/length(spliced_split1@x)

### percentage of non-zero entries
sum(spliced_split1@x>0)/length(spliced_split1@x)
sum(spliced_split2@x>0)/length(spliced_split2@x)

c(spliced_split1@x[which.max(spliced_split1@x)],spliced_split2@x[which.max(spliced_split1@x)])
c(spliced_split1@x[which.max(spliced_split2@x)],spliced_split2@x[which.max(spliced_split2@x)])

## Unpliced splits
### both non-zero
sum(unspliced_split1@x>0 & unspliced_split2@x>0)/length(unspliced_split1@x)

### one of two splits non-zero and the other zero
sum(unspliced_split1@x>0 & unspliced_split2@x==0)/length(unspliced_split1@x)
sum(unspliced_split1@x==0 & unspliced_split2@x>0)/length(unspliced_split1@x)

### percentage of non-zero entries
sum(unspliced_split1@x>0)/length(unspliced_split1@x)
sum(unspliced_split2@x>0)/length(unspliced_split2@x)

c(unspliced_split1@x[which.max(unspliced_split1@x)],unspliced_split2@x[which.max(unspliced_split1@x)])
c(unspliced_split1@x[which.max(unspliced_split2@x)],unspliced_split2@x[which.max(unspliced_split2@x)])



