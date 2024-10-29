rm(list=ls())

library(org.Mm.eg.db)
library(clusterProfiler)

csv_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/csv/kevin/Writeup12/"
df <- read.csv(paste0(csv_folder, "pancreas_gene_likelihood.csv"))
colnames(df)
head(df)
rownames(df) <- df[,"gene_name"]

# one idea: let's just count how many times a gene is highly variable
velocity_var_idx <- grep("X.*velocity_genes", colnames(df))
df_subset <- df[,velocity_var_idx]
num_velocity_var <- apply(df_subset, 1, function(x){length(which(x == "True"))})
names(num_velocity_var) <- df[,"gene_name"]

lik_idx <- grep("lik", colnames(df)) # this is including the total
df_subset <- df[,lik_idx]
median_lik <- apply(df_subset, 1, function(x){stats::median(x, na.rm = TRUE)})
names(median_lik) <- df[,"gene_name"]

# teststat_vec <- num_velocity_var
# teststat_vec <- sort(teststat_vec, decreasing = TRUE)
# keep_genes <- intersect(names(teststat_vec)[teststat_vec != 0],
#                         names(median_lik)[!is.na(median_lik)])
# teststat_vec <- teststat_vec[keep_genes]
# median_lik <- median_lik[keep_genes]
# 
# # break ties in teststat_vec based on the higher likelihood
# uniq_values <- sort(unique(teststat_vec))
# for(val in uniq_values){
#   idx <- which(teststat_vec == val)
#   median_lik_subset <- median_lik[idx]
#   median_lik_subset <- sort(median_lik_subset, decreasing = FALSE)
#   gene_order <- names(median_lik_subset)
#   teststat_vec[gene_order] <- teststat_vec[gene_order] + seq(0, 0.95, length.out = length(gene_order))
# }

teststat_vec <- median_lik
teststat_vec <- teststat_vec[!is.na(teststat_vec)]
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

# https://github.com/YuLab-SMU/clusterProfiler/issues/307
set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 1,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500,            # maximum gene set size
  scoreType = "pos"
)

gse_df <- as.data.frame(gse)

gse_df[which(gse_df$p.adjust <= 0.05), c("Description", "NES", "p.adjust")]
# I don't like these as much... too many neuron things...

write.csv(gse_df, 
          file = paste0(csv_folder, "Writeup12_pancreas_scvelo_GSEA.csv"))


