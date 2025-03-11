rm(list=ls())
library(org.Hs.eg.db)
library(clusterProfiler)

dat <- readxl::read_xlsx("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/Writeup13/1-s2.0-S0092867421009429-mmc2.xlsx",
                         sheet = "B",range = "A2:H1853")
dat <- as.data.frame(dat)

vec <- dat[,"Correlation scaled"]
names(vec) <- dat[,"Gene symbol"]
vec <- sort(vec, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500,            # maximum gene set size
  scoreType = "pos"
)

gse_df <- as.data.frame(gse)

dim(gse_df)
gse_df[1:25,c("ID", "Description")]
gse_df[1:25,c("core_enrichment")]

core_enrichment_genes <- unique(unlist(lapply(1:nrow(gse_df), function(i){
  tmp <- gse_df[i, "core_enrichment"]
  tmp <- strsplit(tmp, split = "/")[[1]]
  tmp
})))

gpc_genes <- dat[dat[,"GPC"] == "Yes","Gene symbol"]

table(gpc_genes %in% core_enrichment_genes)

write.csv(core_enrichment_genes,
          row.names = FALSE,
          "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/Writeup13/Writeup13_blacklist_gsea_gpc_genes.csv")
