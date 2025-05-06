rm(list=ls())

library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(10)

out_folder <- "~/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/"

velovi_woprep <- read.csv(paste0(out_folder, "Writeup15_greenleaf_gene_laplacian_scores_velovi-woprep_transition_neu4.csv"))

# Prepare input for GSEA
teststat_vec <- velovi_woprep$score
names(teststat_vec) <- velovi_woprep$symbol
teststat_vec <- teststat_vec[!is.na(teststat_vec)]
teststat_vec <- sort(teststat_vec, decreasing = TRUE) # Sort in decreasing order

# remove gnes with no symbol
nchar_vec <- sapply(names(teststat_vec), function(x){nchar(x)})
if(any(nchar_vec == 0)) teststat_vec <- teststat_vec[which(nchar_vec > 0)]
  
# Run GSEA
set.seed(10)
gsea_velovi <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # Biological Process
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1
)

gsea_velovi_df <- as.data.frame(gsea_velovi)

dim(gsea_velovi_df)
gsea_velovi_df <- gsea_velovi_df[order(gsea_velovi_df$pvalue, decreasing = FALSE),]
gsea_velovi_df[which(gsea_velovi_df$p.adjust <= 0.05),c("Description", "NES")]
gsea_velovi_df[which(gsea_velovi_df$pvalue <= 0.05),c("Description", "NES", "pvalue")]

######

scvelo <- read.csv(paste0(out_folder, "Writeup15_greenleaf_gene_laplacian_scores_scvelo.csv"))

# Prepare input for GSEA
teststat_vec <- scvelo$score
names(teststat_vec) <- scvelo$symbol
teststat_vec <- teststat_vec[!is.na(teststat_vec)]
teststat_vec <- sort(teststat_vec, decreasing = TRUE) # Sort in decreasing order

# remove gnes with no symbol
nchar_vec <- sapply(names(teststat_vec), function(x){nchar(x)})
if(any(nchar_vec == 0)) teststat_vec <- teststat_vec[which(nchar_vec > 0)]

# Run GSEA
set.seed(10)
gsea_scvelo <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # Biological Process
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  minGSSize = 10,
  maxGSSize = 500
)

gsea_scvelo_df <- as.data.frame(gsea_scvelo)

dim(gsea_scvelo_df)
gsea_scvelo_df[,c("Description", "NES")]


