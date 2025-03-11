rm(list=ls())
library(uwot)

umap_model <- uwot::load_uwot("umapModel_rna_hft_uwot")

names(umap_model)

write.csv(umap_model$embedding, file = "umapModel_rna_hft_uwot.csv")