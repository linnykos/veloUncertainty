rm(list=ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)  

set.seed(10)

plot_folder <- "~/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/fig/kevin/Writeup15/"
out_folder <- "~/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/"
source("~/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/code/kevin/Writeup15_laplacian/plot_signed_logpvalue.R")

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

gsea_df <- gsea_df_list[[3]]
gsea_df$num_core <- sapply(gsea_df$core_enrichment, function(x){
  length(strsplit(x, split = "/")[[1]])
})
idx <- intersect(which(gsea_df$pvalue <= 0.05),
                 which(gsea_df$num_core > 5))
gsea_df_subset <- gsea_df[idx,]
gsea_df_subset[,c("Description", "core_enrichment")]

####################

common_pathways <- intersect(rownames(gsea_df_list[["scvelo"]]),
                             rownames(gsea_df_list[["velovi_woprep"]]))
gsea_velovi_df <- gsea_df_list[["velovi_woprep"]][common_pathways,]
gsea_scvelo_df <- gsea_df_list[["scvelo"]][common_pathways,]

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

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup15_greenleaf_velovo-woprep_scvelo.png"),
                height = 6, 
                width = 6)

####################

go_terms_positive <- c("GO:0051968", "GO:0035249",
                       "GO:0098900", "GO:0043267",
                       "GO:0060070", "GO:0090263", 
                       "GO:0043409", "GO:0046777")

go_terms_negative <- c("GO:0150063", "GO:0048880",
                       "GO:0001654",
                       "GO:0001659", "GO:0001822",
                       "GO:0072001", "GO:0030326")

pvalue_df_positive <- data.frame(sapply(gsea_df_list, function(gsea_df){
  vec <- sapply(go_terms_positive, function(x){
    if(x %in% rownames(gsea_df)){
      return(-log10(gsea_df[x, "pvalue"]))
    } else {
      return(0)
    }
  })
}))
pvalue_df_positive2 <- pvalue_df_positive
pvalue_df_positive2$Description <- gsea_df_list[["scvelo"]][rownames(pvalue_df_positive), "Description"]

pvalue_df_negative <- data.frame(sapply(gsea_df_list, function(gsea_df){
  vec <- sapply(go_terms_negative, function(x){
    if(x %in% rownames(gsea_df)){
      return(-log10(gsea_df[x, "pvalue"]))
    } else {
      return(0)
    }
  })
}))
pvalue_df_negative2 <- pvalue_df_negative
pvalue_df_negative2$Description <- gsea_df_list[["scvelo"]][rownames(pvalue_df_negative), "Description"]


#############

positive_df <- pvalue_df_positive %>%
  rownames_to_column("pathway") %>%                         # keep rownames
  pivot_longer(-pathway, names_to = "method",
               values_to = "log10p") %>%
  mutate(direction = "Positive",
         pathway   = factor(pathway, levels = rownames(pvalue_df_positive))) # keep order

negative_df <- pvalue_df_negative %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "method",
               values_to = "log10p") %>%
  mutate(direction = "Negative",
         pathway   = factor(pathway, levels = rownames(pvalue_df_negative))) # keep order

plot_df <- bind_rows(positive_df, negative_df) %>% 
  mutate(direction = factor(direction,          # ← make it a factor
                            levels = c("Positive", "Negative")))

method_cols <- c(
  sct            = "#1B9E77",  # teal-green  (unchanged)
  scvelo         = "#7570B3",  # indigo      (unchanged)
  utv            = "#666666",  # medium gray with good contrast
  velovi         = "#E6AB02",  # amber       (unchanged)
  velovi_woprep  = "#D55E00"   # saturated red (Okabe-Ito “vermillion”)
)

plot1 <- ggplot(plot_df,
                aes(x = pathway, y = log10p, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  facet_wrap(~ direction, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = method_cols, name = "Method") +
  labs(x = NULL, y = expression(-log[10]~italic(p))) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.y      = element_text(face = "italic"),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.spacing.y  = unit(2, "mm")
  )

ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup15_greenleaf_gsea_barplot.png"),
                height = 6, 
                width = 3)
