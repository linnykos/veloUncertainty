rm(list=ls())

library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(10)

out_folder <- "~/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/"

file_names <- c(
  sct = "Writeup15_greenleaf_gene_laplacian_scores_sct_chung_neu4.csv",
  scvelo = "Writeup15_greenleaf_gene_laplacian_scores_scv_chung_neu4.csv",
  utv = "Writeup15_greenleaf_gene_laplacian_scores_utv_chung_neu4.csv",
  velovi = "Writeup15_greenleaf_gene_laplacian_scores_velovi_chung_neu4.csv",
  velovi_woprep = "Writeup15_greenleaf_gene_laplacian_scores_velovi-woprep_chung_neu4.csv"
)

gsea_df_list <- vector("list", length = length(file_names))
names(gsea_df_list) <- names(file_names)
  
for(kk in 1:length(file_names)){
  method_name <- names(gsea_df_list)[kk]
  file_name <- file_names[kk]
  
  print(paste0("Working on: ", file_name))
  csv_results <- read.csv(paste0(out_folder, file_name))
  
  # Prepare input for GSEA
  teststat_vec <- 1/csv_results$score
  names(teststat_vec) <- csv_results$symbol
  teststat_vec <- teststat_vec[!is.na(teststat_vec)]
  teststat_vec <- sort(teststat_vec, decreasing = TRUE) # Sort in decreasing order
  
  # remove gnes with no symbol
  nchar_vec <- sapply(names(teststat_vec), function(x){nchar(x)})
  if(any(nchar_vec == 0)) teststat_vec <- teststat_vec[which(nchar_vec > 0)]
  
  # Run GSEA
  set.seed(10)
  gsea_results <- clusterProfiler::gseGO(
    teststat_vec,
    ont = "BP", # Biological Process
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    scoreType = "pos"
  )
  
  gsea_df_list[[method_name]] <- as.data.frame(gsea_results)
}

#######

for(i in 1:length(gsea_df_list)){
  print(paste0("Printing: ", names(gsea_df_list)[i]))
  
  gsea_df <- gsea_df_list[[i]]
  gsea_df$num_core <- sapply(gsea_df$core_enrichment, function(x){
    length(strsplit(x, split = "/")[[1]])
  })
  idx <- intersect(which(gsea_df$pvalue <= 0.05),
                   which(gsea_df$num_core > 5))
  print(paste0(length(idx), " number of significant pathways"))
  if(length(idx) > 0){
    print(paste0("Top ", min(length(idx),20), " pathways:"))
    gsea_df_subset <- gsea_df[idx,]
    gsea_df_subset <- gsea_df_subset[order(gsea_df_subset$p.adjust, decreasing = FALSE),]
    print(gsea_df_subset[1:(min(20,length(idx))), "Description"])
  }

  print("====")
}

gsea_df <- gsea_df_list[[5]]
idx <- intersect(which(gsea_df$pvalue <= 0.05),
                 which(gsea_df$num_core > 5))
gsea_df_subset <- gsea_df[idx,]
gsea_df_subset[,c("Description", "core_enrichment")]

####################

common_pathways <- intersect(rownames(gsea_velovi_df),
                             rownames(gsea_scvelo_df))
gsea_velovi_df <- gsea_velovi_df[common_pathways,]
gsea_scvelo_df <- gsea_scvelo_df[common_pathways,]

velovi_idx <- intersect(which(gsea_velovi_df$p.adjust <= 0.05),
                        which(gsea_velovi_df$NES <= 0))
scvelo_idx <- intersect(which(gsea_scvelo_df$p.adjust <= 0.05),
                        which(gsea_scvelo_df$NES <= 0))
length(intersect(velovi_idx, scvelo_idx))
gsea_velovi_df[setdiff(velovi_idx, scvelo_idx), "Description"]
gsea_velovi_df[setdiff(scvelo_idx, velovi_idx), "Description"]

df_combined <- data.frame(
  Gene = common_pathways,
  velovi_logFC = gsea_velovi_df$NES,
  velovi_pval = gsea_velovi_df$pvalue,
  scvelo_logFC = gsea_scvelo_df$NES,
  scvelo_pval = gsea_scvelo_df$pvalue
)

plot1 <- .plot_signed_logpvalue(df_combined,
                                method1 = "velovi",
                                method2 = "scvelo")

plot_folder <- "~/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/fig/kevin/Writeup15/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup15_greenleaf_velovo-woprep_scvelo.png"),
                height = 6, 
                width = 6)

