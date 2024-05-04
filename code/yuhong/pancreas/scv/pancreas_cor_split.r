library(Seurat)
library(SeuratDisk)
library(SeuratObject)

setwd("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split")

split1_seed317 <- LoadH5Seurat("./pancreas_seed317_split1_seurat.h5Seurat")
split2_seed317 <- LoadH5Seurat("./pancreas_seed317_split2_seurat.h5Seurat")
split1_seed320 <- LoadH5Seurat("./pancreas_seed320_split1_seurat.h5Seurat")
split2_seed320 <- LoadH5Seurat("./pancreas_seed320_split2_seurat.h5Seurat")

spliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="spliced"),
                        SeuratObject::LayerData(split2_seed317,assay="spliced"))
unspliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed317,assay="unspliced"))
spliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="spliced"),
                        SeuratObject::LayerData(split2_seed320,assay="spliced"))
unspliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed320,assay="unspliced"))

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

png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/corr_seed317.png"), height = 1000, width = 1800, units = "px", res = 300)
par(mfrow=c(1,2))
plot(cor_spliced_seed317, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Spliced",ylab="corr",cex=.7,cex.main=.8) 
plot(cor_unspliced_seed317, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Unspliced",ylab="corr",cex=.7,cex.main=.8)
graphics.off()

png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/corr_seed320.png"), height = 1000, width = 1800, units = "px", res = 300)
par(mfrow=c(1,2))
plot(cor_spliced_seed320, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Spliced",ylab="corr",cex=.7,cex.main=.8) 
plot(cor_unspliced_seed320, ylim = c(-1,1), pch = 16, col = color_values, 
     main="Unspliced",ylab="corr",cex=.7,cex.main=.8)
graphics.off()



