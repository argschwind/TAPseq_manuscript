## this "snakemake" R-script performs differential gene expression analysis for CROP-seq
## perturbations to discover enhancer - target gene interactions

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# attach required packages and load functions
library(tidyverse)
library(data.table)
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

# get specified method and strategy
method <- snakemake@wildcards$method
strategy <- snakemake@wildcards$strategy

# check if strategy is valid (method is handled by DE functions)
if (!strategy %in% c("perEnh", "perGRNA")) {
  stop("Strategy needs to be either 'perEnh' or 'perGRNA'!", call. = FALSE)
}


# define additional functions ======================================================================

# helper function to convert grna ids to grna target based on information contained in the id
extract_grna_targets <- function(grnas) {
  
  # get indices of screen grnas
  screen_grnas <- grep(grnas, pattern = "^chr.+:\\d+-\\d+.*$", perl = TRUE)
  
  # get non-targeting negative control grnas
  neg_ctrls <- grep(grnas, pattern = "^non-targeting_.+$", perl = TRUE)
  
  # convert grna ids to grna targets
  grnas[screen_grnas] <- sub("_\\d+.+$", "", grnas[screen_grnas])
  grnas[neg_ctrls] <- sub("_\\d+$", "", grnas[neg_ctrls])
  grnas[-c(screen_grnas, neg_ctrls)] <- sub("[_|-].*" , "", grnas[-c(screen_grnas, neg_ctrls)])
  
  # return converted ids
  return(grnas)
  
}


# prepare input ====================================================================================

# load DGE data
dge <- read.table(snakemake@input$dge, header = TRUE, stringsAsFactors = FALSE, row.names = "GENE")

# load perturbation status data
perturb_status <- read.table(snakemake@input$perturb_status, header = TRUE,
                             stringsAsFactors = FALSE)

# rows in dge that contain CROP-seq vector transcripts
vector_rows <- grep(rownames(dge), pattern = snakemake@params$vector_pattern)

# extract expression data on genes
gene_expr <- dge[-vector_rows, ]

# get rows containing negative control perturbations
neg_ctrls <- grep(perturb_status$VECTOR, pattern = "^non-targeting.*$", perl = TRUE)

# remove negative controls, as these are mostly found in control 10x lane
perturb_status <- perturb_status[-neg_ctrls, ]

# remove cells from specified 10x lanes
if (all(snakemake@params$exclude_lanes != "none")) {
  
  # get cells that belong to filtered out lanes
  cells <- colnames(gene_expr)
  cell_lanes <- substr(cells, start = 1, stop = 8)
  cells_filt <- cells[!cell_lanes %in% snakemake@params$exclude_lanes]
  
  # filter perturbation and gene expression data
  perturb_status <- perturb_status[, c("VECTOR", cells_filt)]
  gene_expr <- gene_expr[, cells_filt]
  
}

# convert grna ids to grna targets if specified and count number of cells per perturbation
if (strategy == "perEnh") {
  
  # convert gRNA ids to targeted enhancers
  perturb_status$VECTOR <- extract_grna_targets(perturb_status$VECTOR)
  
  # calculate number of cells per perturbation using data.table (still super slow...)
  count_cells <- function(x) sum(colSums(x) > 0)
  pert_dt <- as.data.table(perturb_status)
  cells_per_pert <- pert_dt[, count_cells(.SD), by = .(VECTOR)]
  cells_per_pert <- data.frame(perturbation = cells_per_pert$VECTOR, cells = cells_per_pert$V1,
                               stringsAsFactors = FALSE)
  
}else{
  
  # simply count cells per gRNA to get cells per perturbation
  cells_per_pert <- rowSums(perturb_status[, -1])
  cells_per_pert <- data.frame(perturbation = perturb_status$VECTOR, cells = cells_per_pert,
                               stringsAsFactors = FALSE)
  
}

# save cells per perturbation to file for later analyses
write.csv(cells_per_pert, file = snakemake@output$cells, row.names = FALSE)


# perform DE tests =================================================================================

# remove perturbations with fewer cells than min_cells
de_perts <- dplyr::filter(cells_per_pert, cells >= snakemake@params$min_cells)
perturb_status <- perturb_status[perturb_status$VECTOR %in% de_perts$perturbation, ]

# remove any cells that do not have any perturbation after filtering
de_cells <- colSums(perturb_status[, -1]) > 0
perturb_status <- perturb_status[, c(TRUE, de_cells)]
gene_expr <- gene_expr[, de_cells]

# perform differential gene expression analysis
output <- test_differential_expression(perturb_status, gene_expr = gene_expr, parallel = parallel,
                                       method = method)

# save output to file
write.csv(output, file = snakemake@output$results, row.names = FALSE)
