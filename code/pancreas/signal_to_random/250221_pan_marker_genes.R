# get all pancreas marker genes
pan_genes <- c(readxl::read_excel('/Users/y2564li/Downloads/proj_scRNA/TableS3.xlsx', sheet='Ductal')$Names,
               readxl::read_excel('/Users/y2564li/Downloads/proj_scRNA/TableS3.xlsx', sheet='Ngn3 low EP')$Names,
               readxl::read_excel('/Users/y2564li/Downloads/proj_scRNA/TableS3.xlsx', sheet='Ngn3 high EP')$Names)
unique(pan_genes)
pan_genes[duplicated(pan_genes)]
length(pan_genes)
length(unique(pan_genes))

df_genes = data.frame('gene_name'=sort(unique(pan_genes)))
write.csv(df_genes, '/Users/y2564li/Downloads/proj_scRNA/pan_marker_gene_name.csv')

