rm(list=ls())

df <- read.csv("~/kzlinlab/data/genelists/housekeeping/HRT_Atlas/Housekeeping_GenesMouse.csv",
               sep = ";")
housekeeping_genes <- sort(unique(df$Gene))

write.table(housekeeping_genes,
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE,
            file = "~/kzlinlab/data/genelists/housekeeping/HRT_Atlas/Housekeeping_GenesMouse_formatted.csv")