rm(list=ls())

library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(10)

out_folder <- "~/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/"

velovi_woprep <- read.csv(paste0(out_folder, "Writeup15_greenleaf_gene_laplacian_scores_velovi-woprep_chung_neu4.csv"))

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
gsea_velovi_df[intersect(which(gsea_velovi_df$pvalue <= 0.05),
                         which(gsea_velovi_df$NES <= 0)),
               c("Description", "NES", "pvalue")]
# GO:0051932                                    synaptic transmission, GABAergic
# GO:0051932 -1.702567 0.024487020

gsea_velovi_df["GO:0051932",]
gsea_velovi_df["GO:1990266",]
GO:0051932


######

scvelo <- read.csv(paste0(out_folder, "Writeup15_greenleaf_gene_laplacian_scores_scv_chung_neu4.csv"))

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
  maxGSSize = 500,
  pvalueCutoff = 1
)

gsea_scvelo_df <- as.data.frame(gsea_scvelo)

dim(gsea_scvelo_df)
gsea_scvelo_df <- gsea_scvelo_df[order(gsea_scvelo_df$pvalue, decreasing = FALSE),]
gsea_scvelo_df[which(gsea_scvelo_df$p.adjust <= 0.05),c("Description", "NES")]
gsea_scvelo_df[intersect(which(gsea_scvelo_df$pvalue <= 0.05),
                         which(gsea_scvelo_df$NES <= 0)),
               c("Description", "NES", "pvalue")]

# GO:1990266                            neutrophil migration -2.051542
# GO:1990266 0.004489973


gsea_scvelo_df["GO:0051932",]
