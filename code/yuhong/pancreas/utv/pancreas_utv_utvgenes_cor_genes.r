library(Seurat)
library(SeuratDisk)
library(SeuratObject)

setwd("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pan_utv_utvgenes")

i_genes_317 <- read.csv('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup9_corr/pan_utv_genes_seed317.csv',header=F)
i_genes_320 <- read.csv('/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup9_corr/pan_utv_genes_seed320.csv',header=F)
i_genes_317 <- sapply(i_genes_317,function(i) i)
i_genes_320 <- sapply(i_genes_320,function(i) i)

gene_names <- rownames(LayerData(split1_seed317))
out_folder_path = "/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/unitvelo_utvgenes/Jun12_cor/"

split1_seed317 <- LoadH5Seurat("./pan_utvgenes_seed317_split1_seurat.h5Seurat")
split2_seed317 <- LoadH5Seurat("./pan_utvgenes_seed317_split2_seurat.h5Seurat")
split1_seed320 <- LoadH5Seurat("./pan_utvgenes_seed320_split1_seurat.h5Seurat")
split2_seed320 <- LoadH5Seurat("./pan_utvgenes_seed320_split2_seurat.h5Seurat")

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

# sum(is.na(cor_spliced_seed317)) = 3
nacor_s317 <- which(is.na(cor_spliced_seed317))
# sum(is.na(cor_unspliced_seed317)) = 22
nacor_u317 <- which(is.na(cor_unspliced_seed317))
# sum(is.na(cor_spliced_seed320)) = 2
nacor_s320 <- which(is.na(cor_spliced_seed320))
# sum(is.na(cor_unspliced_seed320)) = 23
nacor_u320 <- which(is.na(cor_unspliced_seed320))

color_values_317 <- rep("lightblue",2000)
color_values_317[i_genes_317] <- "firebrick"
color_values_320 <- rep("lightblue",2000)
color_values_320[i_genes_320] <- "firebrick"

plot_corr <- function(cor_spliced,cor_unspliced,color_values,figname) {
    png(file=paste0(out_folder_path,figname), height = 1000, width = 2000, units = "px", res = 300)
    par(mfrow=c(1,2))
    plot(cor_spliced, ylim = c(-1,1), pch = 16, col = color_values, main="Spliced",ylab="corr",cex=.7,cex.main=.8) 
    plot(cor_unspliced, ylim = c(-1,1), pch = 16, col = color_values, main="Unspliced",ylab="corr",cex=.7,cex.main=.8)
    graphics.off()
}
plot_corr(cor_spliced=cor_spliced_seed317, cor_unspliced=cor_unspliced_seed317,color_values=color_values_317, figname="corr_seed317.png")
plot_corr(cor_spliced=cor_spliced_seed320, cor_unspliced=cor_unspliced_seed320,color_values=color_values_320, figname="corr_seed320.png")

mean(cor_spliced_seed317,na.rm=T) # 0.03962642
mean(cor_spliced_seed317[i_genes_317]) # 0.03950827
mean(cor_spliced_seed317[setdiff(1:2000,i_genes_317)],na.rm=T) # 0.0396833

plot_corr_hist <- function(i_genes,nacor,cor_spliced,figname) {
    png(file = paste0(out_folder_path,figname), height = 1000, width = 2000, units = "px", res = 300)
    par(mfrow=c(1,2))
    hist(cor_spliced[setdiff(1:2000,c(i_genes,nacor))],col="lightblue",
        main="corr of unselected genes(pan,utv)",cex.main=.8,xlab="correlation")
    abline(v=mean(cor_spliced[setdiff(1:2000,i_genes)],na.rm=T),col="coral",lwd=1.2)
    text(.1,850,label=paste0("mean=",round(mean(cor_spliced[setdiff(1:2000,i_genes)],na.rm=T),4)),col="coral",cex=.5)
    hist(cor_spliced[i_genes],col="firebrick",main="corr of selected genes(pan,utv)",cex.main=.8,xlab="correlation")
    abline(v=mean(cor_spliced[i_genes]),col="coral",lwd=1.2)
    text(.1,130,label=paste0("mean=",round(mean(cor_spliced[i_genes]),4)),col="coral",cex=.5)
    graphics.off()
}
plot_corr_hist(i_genes=i_genes_317,nacor=nacor_s317, cor_spliced=cor_spliced_seed317, figname="corr_seed317_hist_velogenes.png")
plot_corr_hist(i_genes=i_genes_320,nacor=nacor_s320, cor_spliced=cor_spliced_seed320, figname="corr_seed320_hist_velogenes.png")


## plot 24 genes in two spliced count splits with greatest correlations (log10)
intersect(order(abs(cor_spliced_seed317),decreasing=T)[1:25],which(is.na(cor_spliced_seed317))) # integer(0)
intersect(order(abs(cor_spliced_seed317))[1:25],which(is.na(cor_spliced_seed317))) # integer(0)

plot_genes <- function(cor_spliced,spliced,seed) {
    corr_max24 <- order(abs(cor_spliced),decreasing=T)[1:24]
    corr_min24 <- order(abs(cor_spliced),decreasing=F)[1:24]
    for (i_start in as.integer(seq(1,21,length.out=6))) {
        ind = as.integer(seq(i_start,i_start+3,length.out=4)) # 4 genes to be plotted, decreasing corr
        x = list()
        y = list()
        n = length(spliced[[1]][1,])
        for (j in 1:4) {
            i_gene = corr_max24[ind[j]]
            x[[j]] <- log10(spliced[[1]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
            y[[j]] <- log10(spliced[[2]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
        }
        png(file = paste0(out_folder_path,"spliced_seed",seed,"_cor_gene_max",i_start,".png"),
            height = 800, width = 2800, units = "px", res = 300)
        par(mfrow=c(1,4))
        for (j in 1:4) {
            i_gene = corr_max24[ind[j]]
            plot(x[[j]], y[[j]], pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
                main=paste0("gene",i_gene," (",gene_names[i_gene],"), cor=",round(cor_spliced[i_gene],2)),cex=.7,cex.main=.7)
        }
        graphics.off()
    }
    for (i_start in as.integer(seq(1,21,length.out=6))) {
        ind = as.integer(seq(i_start,i_start+3,length.out=4)) # 4 genes to be plotted, decreasing corr
        x = list()
        y = list()
        n = length(spliced[[1]][1,])
        for (j in 1:4) {
            i_gene = corr_min24[ind[j]]
            x[[j]] <- log10(spliced[[1]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
            y[[j]] <- log10(spliced[[2]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
        }
        png(file = paste0(out_folder_path,"spliced_seed",seed,"_cor_gene_min",i_start,".png"),
            height = 800, width = 2800, units = "px", res = 300)
        par(mfrow=c(1,4))
        for (j in 1:4) {
            i_gene = corr_min24[ind[j]]
            plot(x[[j]], y[[j]], pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
                main=paste0("gene",i_gene," (",gene_names[i_gene],"), cor=",round(cor_spliced[i_gene],2)),cex=.7,cex.main=.7)
        }
        graphics.off()
    }
}

plot_genes(cor_spliced=cor_spliced_seed317,spliced=spliced_seed317,seed=317)
plot_genes(cor_spliced=cor_spliced_seed320,spliced=spliced_seed320,seed=320)

