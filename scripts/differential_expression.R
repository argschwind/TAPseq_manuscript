## this "snakemake" R-script performs differential gene expression analysis for CROP-seq
## perturbations to discover enhancer - target gene interactions

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required functions
source("scripts/SingleCellExperiment_tapseq.R")
source("scripts/differential_expression_fun.R")

# register parallel backend if specified (if more than 1 thread provided)
parallel <- FALSE
if (snakemake@threads > 1) {
  library(BiocParallel)
  register(MulticoreParam(workers = snakemake@threads))
  parallel <- TRUE
}

# set seed for reproducible results
if (snakemake@params$seed > 0) set.seed(snakemake@params$seed)

# prepare input data ===============================================================================

# load DGE data
dge <- read.table(snakemake@input$dge, header = TRUE, stringsAsFactors = FALSE, row.names = "GENE")

# load perturbation status data
perturb_status <- read.table(snakemake@input$perturb_status, header = TRUE,
                             stringsAsFactors = FALSE, row.names = 1)

# remove negative controls, as these are not supposed to be tested
neg_ctrls <- grep(rownames(perturb_status), pattern = "^non-targeting.*$", perl = TRUE)
perturb_status <- perturb_status[-neg_ctrls, ]

# count cells per perturbation and save to file for later analyses
cells_per_pert <- rowSums(perturb_status > 0)
cells_per_pert <- data.frame(perturbation = rownames(perturb_status), cells = cells_per_pert)
write.csv(cells_per_pert, file = snakemake@output$ncells, row.names = FALSE)

# create tap-seq SingleCellExperiment object
sce <- create_tapseq_sce(dge = dge, perturb_status = perturb_status,
                         vector_pattern = snakemake@params$vector_pattern)

# remove cells from specified 10x lanes
if (all(snakemake@params$exclude_lanes != "none")) {
  
  # get cells that do not belong to filtered out lanes
  cells <- colnames(sce)
  cell_lanes <- substr(cells, start = 1, stop = 8)
  cells_filt <- cells[!cell_lanes %in% snakemake@params$exclude_lanes]
  
  # only retain data on filtered cells
  sce <- sce[, cells_filt]
  
}

# perform PCA on gene expression data, extract PC loadings, and add as colData to use as covariates
pca_genex <- prcomp(t(assay(sce)), scale. = TRUE)
pcs <- seq_len(as.integer(snakemake@wildcards$npcs))
pc_loadings <- DataFrame(pca_genex$x[, pcs])
colData(sce) <- pc_loadings

# perform DE tests =================================================================================

# perform differential gene expression analysis
output <- test_differential_expression(sce, min_cells = snakemake@params$min_cells,
                                       method = "MAST", parallel = parallel)

# save output
write.csv(output, file = snakemake@output$results, row.names = FALSE)
