setwd("/Users/wakeup/Downloads/UW_23-25/proj_RNA/data/erythroid")

split1_seed317 <- LoadH5Seurat("./split_data/scvelo_erythroid_split1_seurat_seed317.h5Seurat")
split2_seed317 <- LoadH5Seurat("./split_data/scvelo_erythroid_split2_seurat_seed317.h5Seurat")
spliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="spliced"),
                        SeuratObject::LayerData(split2_seed317,assay="spliced"))
unspliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed317,assay="unspliced"))

cor_spliced_seed317 <- sapply(1:nrow(spliced_seed317[[1]]),function(i){ 
  cor( log10(spliced_seed317[[1]][i,]+1), log10(spliced_seed317[[2]][i,]+1) ) } )

cor_unspliced_seed317 <- sapply(1:nrow(unspliced_seed317[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(unspliced_seed317[[1]][i,]+1), log10(unspliced_seed317[[2]][i,]+1) ) } )

#######################
# random code by kevin
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
plot(cor_spliced_seed317, ylim = c(-1,1), pch = 16, col = color_values) # gray represents a small fraction of zeros
plot(cor_unspliced_seed317, ylim = c(-1,1), pch = 16, col = color_values)

i <- which.min(cor_spliced_seed317)
hist(total[i,])

# overlaying the negative binomial dist
tmp <- total[i,]
res <- MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(tmp)))
overdisp <- res$theta
x_vec <- round(seq(0, 120, length.out = 100))
y_vec <- stats::dnbinom(x_vec, size = overdisp, mu = mean(tmp))

hist(total[i,], main = "Overlay NB")
max_val <- 3500
y_vec <- y_vec * max_val/max(y_vec)
points(x_vec, y_vec, col = "coral", pch = 16)

# overlaying the poisson dist
tmp <- total[i,]
x_vec <- round(seq(0, 120, length.out = 100))
y_vec <- stats::dpois(x_vec, lambda = mean(tmp))

hist(total[i,], main = "Overlay Poisson")
max_val <- 3500
y_vec <- y_vec * max_val/max(y_vec)
points(x_vec, y_vec, col = "coral", pch = 16)

i = which.min(abs(cor_spliced_seed317)) # 1463
cor_spliced_seed317[i]
plot(log10(spliced_seed317[[1]][i,]+1),log10(spliced_seed317[[2]][i,]+1),pch=16,asp=T,col=rgb(0.5,0.5,0.5,0.5))
## jittered plots
x = log10(spliced_seed317[[1]][i,]+1)
y = log10(spliced_seed317[[2]][i,]+1)
n = length(x)
x = x + runif(n, min = -0.1, max = 0.1)
y = y + runif(n, min = -0.1, max = 0.1)
plot(x, y,pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1))


#####################
## multiple genes (spliced counts)
head(order(abs(cor_spliced_seed317),decreasing=T)) # 70  630 1649 1010 1552 1574

png(file = "./spliced_seed317_cor.png", height = 1000, width = 1100, units = "px", res = 300)
par(mar = c(4, 4, 1, 1))
plot(cor_spliced_seed317, ylim = c(-1,1), pch = 16, col = color_values)
graphics.off()

## plot counts in two splits for 6 genes with the most strong negative correlations
png(file = paste0("./spliced_seed317_negcor_genes.png"),
    height = 1500, width = 2200, units = "px", res = 300)
par(mfrow=c(2,3))
for (i in head(order(abs(cor_spliced_seed317),decreasing=T)) ) {
  x = log10(spliced_seed317[[1]][i,]+1)
  y = log10(spliced_seed317[[2]][i,]+1)
  n = length(x)
  x = x + runif(n, min = -0.1, max = 0.1)
  y = y + runif(n, min = -0.1, max = 0.1)
  plot(x, y,pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
       main=paste0("gene",i," log10, jittered"),cex=.7)
}
graphics.off()

# overlaying the negative binomial dist
png(file = paste0("./spliced_seed317_overlayNB.png"),
    height = 1500, width = 2200, units = "px", res = 300)
par(mfrow=c(2,3))
for (i in head(order(abs(cor_spliced_seed317),decreasing=T)) ) {
  tmp <- total[i,]
  res <- MASS::glm.nb(x ~ 1, data=data.frame(x=as.numeric(tmp)))
  overdisp <- res$theta
  x_vec <- round(seq(0, 120, length.out = 100))
  y_vec <- stats::dnbinom(x_vec, size = overdisp, mu = mean(tmp))
  hist(total[i,], col="lightblue",main=paste0("total counts - gene",i),xlab="")
  max_val <- 3000
  y_vec <- y_vec * max_val/max(y_vec)
  lines(x_vec, y_vec, col = "coral",lwd=1.5)
}
graphics.off()


### unspliced counts
png(file = "./unspliced_seed317_cor.png", height = 1000, width = 1100, units = "px", res = 300)
par(mar = c(4, 4, 1, 1))
plot(cor_unspliced_seed317, ylim = c(-1,1), pch = 16, col = color_values)
graphics.off()

png(file = paste0("./unspliced_seed317_negcor_genes.png"),
    height = 1500, width = 2200, units = "px", res = 300)
par(mfrow=c(2,3))
for (i in head(order(abs(cor_unspliced_seed317),decreasing=T)) ) {
  x = log10(unspliced_seed317[[1]][i,]+1)
  y = log10(unspliced_seed317[[2]][i,]+1)
  n = length(x)
  x = x + runif(n, min = -0.1, max = 0.1)
  y = y + runif(n, min = -0.1, max = 0.1)
  plot(x, y,pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
       main=paste0("gene",i," log10, jittered"),cex=.7)
}
graphics.off()





