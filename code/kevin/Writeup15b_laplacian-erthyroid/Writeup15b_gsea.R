rm(list=ls())

library(clusterProfiler)
library(org.Mm.eg.db)

set.seed(10)

out_folder <- "~/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15b/"

file_names <- c(
  scvelo = "Writeup15b_erythroid_gene_laplacian_scores_utv_chung.csv"
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
  names(teststat_vec) <- csv_results$gene
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
    OrgDb = "org.Mm.eg.db",
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
  idx <- which(gsea_df$p.adjust <= 0.05)
  print(paste0(length(idx), " number of significant pathways"))
  if(length(idx) > 0){
    print(paste0("Top ", min(length(idx),20), " pathways:"))
    gsea_df_subset <- gsea_df[idx,]
    gsea_df_subset <- gsea_df_subset[order(gsea_df_subset$p.adjust, decreasing = FALSE),]
    print(gsea_df_subset[1:(min(20,length(idx))), "Description"])
  }

  print("====")
}
