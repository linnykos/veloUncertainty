rm(list=ls())

library(Seurat)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/Writeup14/pijunsala-erthyroid.RData")
seurat_obj
head(seurat_obj@meta.data)

kucinski_df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/csv/kevin/Writeup14/Writeup14_kucinski_erthyroid_genes.csv")
head(kucinski_df)

table(kucinski_df$X %in% SeuratObject::Features(seurat_obj))
kucinski_genes <- kucinski_df$X[kucinski_df$X %in% SeuratObject::Features(seurat_obj)]

Seurat::DimPlot(seurat_obj,
                group.by = "celltype",
                cols = seurat_obj@misc$celltype_colors)


#########################

library(princurve)

# do a simple principle curve analysis
princurve_res <- princurve::principal_curve(pca_mat[,1:30])
seurat_obj$pseudotime <- princurve_res$lambda


scCustomize::FeaturePlot_scCustom(seurat_obj, 
                                  features = "pseudotime", 
                                  reduction = "umap")



mat <- as.matrix(SeuratObject::LayerData(seurat_obj,
                                         layer = "data",
                                         assay = "RNA"))
corr_mat <- sapply(1:nrow(mat), function(j){
  if(j %% floor(nrow(mat)/10) == 0) cat('*')
  
  tmp <- stats::cor.test(x = seurat_obj$pseudotime,
                         y = mat[j,])
  c(value = tmp$estimate, pvalue = tmp$p.value)
})
corr_mat <- t(corr_mat)
corr_mat <- as.data.frame(corr_mat)
rownames(corr_mat) <- rownames(mat)

###########

corr_values <- corr_mat[kucinski_genes, "value.cor"]  # Select relevant genes
corr_values <- round(corr_values, 2)  # Round to 2 significant digits
names(corr_values) <- kucinski_genes


# Define output PDF file
pdf_file <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/fig/kevin/Writeup14/Writeup14_pijunsala_erthyroid_genes.pdf"  # Change to your desired path
pdf(pdf_file, width = 10, height = 10, onefile = TRUE)  # Adjust size as needed

num_plots <- 3
# Loop through genes in groups of 25 per page
for (i in 1:floor(length(kucinski_genes)/num_plots^2)) {
  print(i)
  gene_subset <- kucinski_genes[((i-1)*num_plots^2+1):min(i*num_plots^2, length(kucinski_genes))]  # Get up to 25 genes
  
  # Generate FeaturePlots
  plots <- lapply(gene_subset, function(gene) {
    # Get correlation value for the gene
    cor_val <- ifelse(!is.na(corr_values[gene]), paste0(" (r = ", corr_values[gene], ")"), "")
    
    # Create the FeaturePlot
    p <- scCustomize::FeaturePlot_scCustom(seurat_obj, features = gene, reduction = "umap", raster = TRUE) +
      ggtitle(paste0(gene, cor_val)) +  # Title with gene and correlation
      theme(axis.text = element_blank(),  # Remove x and y axis labels
            axis.ticks = element_blank(),
            axis.title = element_blank())
    return(p)
  })
  
  print(plot_grid(plotlist = plots, ncol = num_plots, nrow = num_plots))  # 5x5 grid per page
}

dev.off()  # Close PDF file


