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

total = spliced_seed317[[1]] + spliced_seed317[[2]]
p <- nrow(total)
zero_fraction = sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  length(which(total[j,] == 0))/ncol(total)
}) ## fraction of cells with zero counts for each gene
break_values <- seq(0, 1, length.out = 20)
color_palette <- grDevices::colorRampPalette(c("lightgray", "coral"))(20)
color_values <- sapply(zero_fraction, function(val){
  color_palette[which.min(abs(val - break_values))]
})

cor_spliced_seed317 <- sapply(1:nrow(spliced_seed317[[1]]),function(i){ 
  cor( log10(spliced_seed317[[1]][i,]+1), log10(spliced_seed317[[2]][i,]+1) ) } )

cor_unspliced_seed317 <- sapply(1:nrow(unspliced_seed317[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(unspliced_seed317[[1]][i,]+1), log10(unspliced_seed317[[2]][i,]+1) ) } )



#####################
## multiple genes (spliced counts)

## plot counts in two splits for 6 genes with the most strong negative correlations
for (i in order(abs(cor_spliced_seed317),decreasing=T)[1:2] ) {
  x = log10(spliced_seed317[[1]][i,]+1)
  y = log10(spliced_seed317[[2]][i,]+1)
  n = length(x)
  x = x + runif(n, min = -0.1, max = 0.1)
  y = y + runif(n, min = -0.1, max = 0.1)
  png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/spliced_seed317_negcor_gene",i,".png"),
    height = 1000, width = 1800, units = "px", res = 300)
  plot(x, y,pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
       main=paste0("gene",i," log10, jittered"),cex=.7)
}
graphics.off()

# overlaying the negative binomial dist
for (i in order(abs(cor_spliced_seed317),decreasing=T)[1:2] ) {
  tmp <- total[i,]
  res <- MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(tmp)))
  overdisp <- res$theta
  x_vec <- round(seq(0, 120, length.out = 100))
  y_vec <- stats::dnbinom(x_vec, size = overdisp, mu = mean(tmp))
  png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/spliced_seed317_overlayNB_gene",i".png"),height = 1000, width = 1800, units = "px", res = 300)
  hist(total[i,], col="lightblue",main=paste0("total counts - gene",i),xlab="")
  max_val <- 3000
  y_vec <- y_vec * max_val/max(y_vec)
  lines(x_vec, y_vec, col = "coral",lwd=1.5)
}
graphics.off()


