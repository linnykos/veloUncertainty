rm(list=ls())
library(Seurat)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/Writeup14/kucinski_seurat.RData")

Seurat::DimPlot(seurat_obj,
                reduction = "umap",
                group.by = "anno_man", 
                label = TRUE)

Seurat::FeaturePlot(seurat_obj,
                    reduction = "umap",
                    features = "dpt_pseudotime")

seurat_obj <- subset(seurat_obj, anno_man == "Ery")

mat <- SeuratObject::LayerData(seurat_obj,
                               layer = "scale.data",
                               assay = "RNA")
sd_vec <- matrixStats::rowSds(mat)
mat <- mat[sd_vec >= 1e-3,]
pseudotime <- seurat_obj$dpt_pseudotime

corr_mat <- sapply(1:nrow(mat), function(j){
  if(j %% floor(nrow(mat)/10) == 0) cat('*')
  
  tmp <- stats::cor.test(x = pseudotime,
                         y = mat[j,])
  c(value = tmp$estimate, pvalue = tmp$p.value)
})
corr_mat <- t(corr_mat)
corr_mat <- as.data.frame(corr_mat)
rownames(corr_mat) <- rownames(mat)

par(mfrow = c(1,2))
hist(-log10(corr_mat$pvalue))
plot(sort(corr_mat$value.cor))

#############

library(clusterProfiler)
library(org.Mm.eg.db)

teststat_vec <- corr_mat$value.cor
names(teststat_vec) <- rownames(corr_mat)
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500,            # maximum gene set size
  # scoreType = "pos",
  eps = 0
)

gse_df <- as.data.frame(gse)
dim(gse_df)
gse_df[,c("ID", "Description")]

gse_df[c("GO:0030097", "GO:1903706"),]

#######

set.seed(10)
gse_all <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 1,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500,            # maximum gene set size
  # scoreType = "pos",
  eps = 0
)

gse_all_df <- as.data.frame(gse_all)
dim(gse_all_df)

which(c("GO:0030218", "GO:0048821", "GO:0060319") %in% rownames(gse_all_df))
gse_all_df["GO:0030218",]

#######

erythrocyte_differentiation <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/csv/kevin/Writeup14/GO_term_summary_GO-0030218_erythrocyte-differentiation.txt", 
                                        sep = "\t", row.names = NULL)$MGI.Gene.Marker.ID
erythrocyte_development <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/csv/kevin/Writeup14/GO_term_summary_GO-0048821_erythrocyte-development.txt", 
                                    sep = "\t", row.names = NULL)$MGI.Gene.Marker.ID
primitive_erythrocyte_differentation <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/csv/kevin/Writeup14/GO_term_summary_GO-0060319_primitive-erythrocyte-differentation.txt", 
                                                 sep = "\t", row.names = NULL)$MGI.Gene.Marker.ID

# Define your custom GO pathways
my_go_list <- list(
  "erythrocyte_differentiation" = unique(erythrocyte_differentiation),
  "erythrocyte_development" = unique(erythrocyte_development),
  "primitive_erythrocyte_differentation" = unique(primitive_erythrocyte_differentation)
)

my_go_list <- sapply(my_go_list, function(x){
  x[x %in% names(teststat_vec)]
})

TERM2GENE <- do.call(rbind, lapply(names(my_go_list), function(term) {
  data.frame(term, gene = my_go_list[[term]], stringsAsFactors = FALSE)
}))

# Run GSEA with your custom pathways
gse_custom <- clusterProfiler::GSEA(
  teststat_vec,
  TERM2GENE = TERM2GENE,
  pvalueCutoff = 1,
  minGSSize = 0,            # minimum gene set size
  maxGSSize = 500, 
)
gse_custom_df <- as.data.frame(gse_custom)
gse_custom_df

#############################

# collect all the genes
gene_list1 <- apply(gse_df[c("GO:0030097", "GO:1903706"),], 1, function(x){
  tmp <- x["core_enrichment"]
  strsplit(tmp, split = "/")[[1]]
})
gene_list2 <- strsplit(gse_all_df["GO:0030218","core_enrichment"], split = "/")
gene_list3 <- apply(gse_custom_df, 1, function(x){
  tmp <- x["core_enrichment"]
  strsplit(tmp, split = "/")[[1]]
})

gene_vec <- unique(c(unlist(gene_list1), unlist(gene_list2), unlist(gene_list3)))

hist(-log10(corr_mat$pvalue))
rug(-log10(corr_mat[gene_vec,"pvalue"]), lwd = 2, col = 2)

###################

library(ggplot2)
library(cowplot)
library(gridExtra)

# Extract the correlation values for each gene
corr_values <- corr_mat[gene_vec, "value.cor"]  # Select relevant genes
corr_values <- round(corr_values, 2)  # Round to 2 significant digits
names(corr_values) <- gene_vec

# Define output PDF file
pdf_file <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/fig/kevin/Writeup14/Writeup14_kucinski_erthyroid_genes.pdf"  # Change to your desired path
pdf(pdf_file, width = 10, height = 10, onefile = TRUE)  # Adjust size as needed

num_plots <- 3
# Loop through genes in groups of 25 per page
for (i in 1:floor(length(gene_vec)/num_plots^2)) {
  print(i)
  gene_subset <- gene_vec[((i-1)*num_plots^2+1):min(i*num_plots^2, length(gene_vec))]  # Get up to 25 genes
  
  # Generate FeaturePlots
  plots <- lapply(gene_subset, function(gene) {
    # Get correlation value for the gene
    cor_val <- ifelse(!is.na(corr_values[gene]), paste0(" (r = ", corr_values[gene], ")"), "")
    
    # Create the FeaturePlot
    p <- FeaturePlot(seurat_obj, features = gene, reduction = "umap", raster = TRUE) +
      ggtitle(paste0(gene, cor_val)) +  # Title with gene and correlation
      theme(axis.text = element_blank(),  # Remove x and y axis labels
            axis.ticks = element_blank(),
            axis.title = element_blank())
    return(p)
  })
  
  print(plot_grid(plotlist = plots, ncol = num_plots, nrow = num_plots))  # 5x5 grid per page
}

dev.off()  # Close PDF file

######

corr_mat_subset <- corr_mat[gene_vec,]
corr_mat_subset <- corr_mat_subset[order(abs(corr_mat_subset$value.cor), decreasing = TRUE),]

write.csv(corr_mat_subset, 
          file = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/csv/kevin/Writeup14/Writeup14_kucinski_erthyroid_genes.csv")
