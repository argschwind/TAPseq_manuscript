## script to filter mouse cell type mix umi observations for genes that are expressed in all cells

# load input files
umi_obs <- read.table(snakemake@input$umi_obs, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_genes <- read.csv(snakemake@input$target_genes, stringsAsFactors = FALSE)

# get genes that are expressed in all cell types
genes_all_cells <- target_genes[target_genes$marker == "all", "gene"]

# get umi observations for these genes
umi_obs_out <- umi_obs[umi_obs$Gene %in% genes_all_cells, ]

# save to output file
colnames(umi_obs_out)[1] <- "Cell Barcode"
write.table(umi_obs_out, file = snakemake@output[[1]], sep = "\t", row.names = FALSE, quote = FALSE)
