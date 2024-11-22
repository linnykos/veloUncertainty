rm(list=ls())

library(tidyverse)

csv_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/csv/kevin/Writeup12/"

gsea_velovi <- read.csv(paste0(csv_folder, "Writeup12_pancreas_velovi_GSEA.csv"))
rownames(gsea_velovi) <- gsea_velovi$X
gsea_unitvelo <- read.csv(paste0(csv_folder, "Writeup12_pancreas_unitvelo_GSEA.csv"))
rownames(gsea_unitvelo) <- gsea_unitvelo$X
gsea_scvelo <- read.csv(paste0(csv_folder, "Writeup12_pancreas_scvelo_GSEA.csv"))
rownames(gsea_scvelo) <- gsea_scvelo$X

# find if there's any intersection in terms
threshold <- 0.05
go_velovi_sig <- gsea_velovi$ID[gsea_velovi$p.adjust <= threshold]
go_unitvelo_sig <- gsea_unitvelo$ID[gsea_unitvelo$p.adjust <= threshold]
go_scvelo_sig <- gsea_scvelo$ID[gsea_scvelo$p.adjust <= threshold]

print(length(go_velovi_sig))
print(length(go_unitvelo_sig))
print(length(go_scvelo_sig))

gsea_velovi[go_velovi_sig, "Description"]
gsea_unitvelo[go_unitvelo_sig, "Description"]

######

threshold <- 0.1
go_velovi_sig <- gsea_velovi$ID[gsea_velovi$p.adjust <= threshold]
go_unitvelo_sig <- gsea_unitvelo$ID[gsea_unitvelo$p.adjust <= threshold]
go_scvelo_sig <- gsea_scvelo$ID[gsea_scvelo$p.adjust <= threshold]

print(length(go_velovi_sig))
print(length(go_unitvelo_sig))
print(length(go_scvelo_sig))

table(table(c(go_velovi_sig, go_unitvelo_sig, go_scvelo_sig)))
go_intersect <- intersect(
  intersect(go_velovi_sig, go_unitvelo_sig),
  go_scvelo_sig
)

# find what go terms these are
gsea_velovi[go_intersect, "Description"]

gsea_velovi[go_intersect, c("NES", "p.adjust")]
gsea_unitvelo[go_intersect, c("NES", "p.adjust")]
gsea_scvelo[go_intersect, c("NES", "p.adjust")]

gsea_velovi[1,c("Description", "p.adjust")]
gsea_unitvelo[which(gsea_unitvelo$p.adjust == min(gsea_unitvelo$p.adjust)),
              c("Description", "p.adjust")]

########################

df <- data.frame(id = gsea_velovi[go_intersect, "ID"],
                 Description = gsea_velovi[go_intersect, "Description"],
                 scvelo = -log10(gsea_scvelo[go_intersect, "p.adjust"]),
                 unitvelo = -log10(gsea_unitvelo[go_intersect, "p.adjust"]),
                 velovi = -log10(gsea_velovi[go_intersect, "p.adjust"]))

library(tidyr)

df_long <- df %>%
  pivot_longer(
    cols = c(scvelo, unitvelo, velovi),
    names_to = "method",
    values_to = "logpvalue"
  )

# colors
#E69F00 - Golden orange
#56B4E9 - Sky blue
#009E73 - Bluish green
#F0E442 - Yellow
#CC79A7 - Pink

df_long$method <- factor(df_long$method, levels = c("scvelo", "unitvelo", "velovi"))
plot1 <- ggplot(df_long, aes(x = Description, y = logpvalue, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +  # "dodge" to compare side-by-side bars
  labs(x = "Pathways", y = "-log10(p-adjust)", 
       title = "-log10(p-adjust) of Gene Pathways") + scale_fill_manual(values = c(
         "scvelo" = "#E69F00",  
         "unitvelo" = "#56B4E9",  
         "velovi" = "#009E73"   
       )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/fig/kevin/Writeup12/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup12_gsea_horizontal_barplot.png"),
                height = 1200, width = 3000, units = "px")

##################

.process_string <- function(vec){
  vec <- sapply(vec, function(x){
    strsplit(x, split = "/")[[1]]
  })
  unique(unlist(vec))
}

# find the genes in any of these sets:
scvelo_genes <- .process_string(gsea_scvelo[go_intersect, "core_enrichment"])
unitvelo_genes <- .process_string(gsea_unitvelo[go_intersect, "core_enrichment"])
velovi_genes <- .process_string(gsea_velovi[go_intersect, "core_enrichment"])

tab_vec <- table(c(scvelo_genes, unitvelo_genes, velovi_genes))
enriched_genes <- unique(c(scvelo_genes, unitvelo_genes, velovi_genes))

# start with scvelo
df <- read.csv(paste0(csv_folder, "pancreas_gene_likelihood.csv"))
velocity_var_idx <- grep("X.*velocity_genes", colnames(df))
df_subset <- df[,velocity_var_idx]
num_velocity_var <- apply(df_subset, 1, function(x){length(which(x == "True"))})
names(num_velocity_var) <- df[,"gene_name"]

lik_idx <- grep("X.*lik", colnames(df)) # this is including the total
df_subset <- df[,lik_idx]
median_lik <- apply(df_subset, 1, function(x){stats::median(x, na.rm = TRUE)})
names(median_lik) <- df[,"gene_name"]

teststat_vec <- num_velocity_var
teststat_vec <- sort(teststat_vec, decreasing = TRUE)
keep_genes <- intersect(names(teststat_vec)[teststat_vec != 0],
                        names(median_lik)[!is.na(median_lik)])
teststat_vec <- teststat_vec[keep_genes]
median_lik <- median_lik[keep_genes]

df <- data.frame(gene = keep_genes,
                 count = teststat_vec,
                 median_lik = median_lik,
                 marked = keep_genes %in% enriched_genes)
df <- df %>%
  arrange(desc(count), desc(median_lik))
# Select the top 30 genes based on the sorted data frame
df_top <- df %>% slice(1:50) %>%
  mutate(gene = paste0(gene, " (", count, ")"))

# Create the horizontal bar plot
plot1 <- ggplot(df_top, aes(x = median_lik, y = reorder(gene, median_lik), fill = marked)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "darkgray")) +
  labs(x = "Median likelihood", y = "Gene", 
       title = "Top 50 Genes by Median Likelihood") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = ifelse(rev(df_top$marked), "red", "darkgray")),
    legend.position = "none"
  )

plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/fig/kevin/Writeup12/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup12_scvelo-genes.png"),
                height = 1800, width = 800, units = "px")

#################

# velovi
df <- read.csv(paste0(csv_folder, "pancreas_velovi_gene_velocityR2.csv"))
velocity_var_idx <- grep("X.*velocity_genes", colnames(df))
df_subset <- df[,velocity_var_idx]
num_velocity_var <- apply(df_subset, 1, function(x){length(which(x == "True"))})
names(num_velocity_var) <- df[,"gene_name"]

lik_idx <- grep("X.*velocity_r2", colnames(df)) # this is including the total
df_subset <- df[,lik_idx]
median_lik <- apply(df_subset, 1, function(x){stats::median(x, na.rm = TRUE)})
names(median_lik) <- df[,"gene_name"]

teststat_vec <- num_velocity_var
teststat_vec <- sort(teststat_vec, decreasing = TRUE)
keep_genes <- intersect(names(teststat_vec)[teststat_vec != 0],
                        names(median_lik)[!is.na(median_lik)])
teststat_vec <- teststat_vec[keep_genes]
median_lik <- median_lik[keep_genes]

df <- data.frame(gene = keep_genes,
                 count = teststat_vec,
                 median_lik = median_lik,
                 marked = keep_genes %in% enriched_genes)
df <- df %>%
  arrange(desc(count), desc(median_lik))
# Select the top 30 genes based on the sorted data frame
df_top <- df %>% slice(1:50) %>%
  mutate(gene = paste0(gene, " (", count, ")"))

# Create the horizontal bar plot
plot1 <- ggplot(df_top, aes(x = median_lik, y = reorder(gene, median_lik), fill = marked)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "darkgray")) +
  labs(x = "Median likelihood", y = "Gene", 
       title = "Top 50 Genes by Median Likelihood") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = ifelse(rev(df_top$marked), "red", "darkgray")),
    legend.position = "none"
  )

plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/fig/kevin/Writeup12/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup12_velovi-genes.png"),
                height = 1800, width = 800, units = "px")

###########

# unitvelo
df <- read.csv(paste0(csv_folder, "pancreas_utv_gene_fitloss.csv"))
velocity_var_idx <- grep("X.*velocity_genes", colnames(df))
df_subset <- df[,velocity_var_idx]
num_velocity_var <- apply(df_subset, 1, function(x){length(which(x == "True"))})
names(num_velocity_var) <- df[,"gene_name"]

lik_idx <- grep("X.*fit_loss", colnames(df)) # this is including the total
df_subset <- df[,lik_idx]
median_lik <- apply(df_subset, 1, function(x){stats::median(x, na.rm = TRUE)})
names(median_lik) <- df[,"gene_name"]

teststat_vec <- num_velocity_var
teststat_vec <- sort(teststat_vec, decreasing = TRUE)
keep_genes <- intersect(names(teststat_vec)[teststat_vec != 0],
                        names(median_lik)[!is.na(median_lik)])
teststat_vec <- teststat_vec[keep_genes]
median_lik <- median_lik[keep_genes]

df <- data.frame(gene = keep_genes,
                 count = teststat_vec,
                 median_lik = median_lik,
                 marked = keep_genes %in% enriched_genes)
df <- df %>%
  arrange(desc(count), desc(median_lik))
# Select the top 30 genes based on the sorted data frame
df_top <- df %>% slice(1:50) %>%
  mutate(gene = paste0(gene, " (", count, ")"))

# Create the horizontal bar plot
plot1 <- ggplot(df_top, aes(x = median_lik, y = reorder(gene, median_lik), fill = marked)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "darkgray")) +
  labs(x = "Median likelihood", y = "Gene", 
       title = "Top 50 Genes by Median Likelihood") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = ifelse(rev(df_top$marked), "red", "darkgray")),
    legend.position = "none"
  )

plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/git/veloUncertainty/fig/kevin/Writeup12/"
ggplot2::ggsave(plot1, 
                filename = paste0(plot_folder, "Writeup12_unitvelo-genes.png"),
                height = 1800, width = 800, units = "px")
