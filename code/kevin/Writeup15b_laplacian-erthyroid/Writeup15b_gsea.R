rm(list=ls())

library(clusterProfiler)
library(org.Mm.eg.db)

set.seed(10)

plot_folder <- "~/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/fig/kevin/Writeup15b/"
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
  
  gsea_df_list[[method_name]] <- gsea_results
}

#######

for(i in 1:length(gsea_df_list)){
  print(paste0("Printing: ", names(gsea_df_list)[i]))
  
  gsea_df <- data.frame(gsea_df_list[[i]])
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

########

gsea_results <- gsea_df_list[[1]]

res_tbl <- gsea_results@result %>%                     # pull out the result slot
  mutate(coreSize = str_count(core_enrichment, "/") + 1) %>%   # genes per set
  filter(pvalue <= 0.05, coreSize > 5) %>%             # your two rules
  arrange(pvalue) %>%                                  # optional: sort
  slice_head(n = 20)                                   # keep top-20 after filter

gsea_top20 <- gsea_results          # shallow copy
gsea_top20@result <- res_tbl        # replace the internal table

plot1 <- clusterProfiler::dotplot(gsea_top20,                # still a gseaResult object
                                  showCategory = 20,         # number of rows now in @result
                                  orderBy      = "pvalue",
                                  color = "pvalue") +
  ggtitle("Top 20 filtered GSEA pathways\n(Erythroid, UnitVelo)")

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup15b_erythroid_gsea_dotplot.png"),
                height = 12, 
                width = 8)


