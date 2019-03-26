## functions to perform differential gene expression analysis for enhancer TAP-seq ##

library(Seurat)

#' Censored mean normalization
#' 
#' Normalize DGE data based on censored mean, which excludes genes highly expressed genes when
#' calculating size factors.
#' 
#' @param dge Matrix or data.frame containing digital gene expression data for all cells and genes
#'   to be tested. Rows are genes and columns are cells, with row names being gene ids and column
#'   names cell ids.
#' @param percentile Genes of each cell with expression from that percentile upwards will be
#'   excluded to calculate size factors.
#' @param norm_counts (logical) Return normalized counts, meaning normalized data is rounded to the
#'   nearest integer to resemble a raw count matrix required by certain DE methods such as DEsingle.
#'   Default = FALSE.
#' @param scale_fun Scaling function to be applied to censored mean normalized data. Ignored if
#'   norm_counts = TRUE.
normalize_cens_mean <- function(dge, percentile = 0.9, norm_counts = FALSE, scale_fun = log1p) {
  
  # function to calculate censored size factor for one cell
  calculate_sf <- function(x) {
    sum(x[x <= quantile(x, probs = percentile)]) + 1
  }
  
  # calculate size factors for normalization
  size_factors <- apply(dge, MARGIN = 2, FUN = calculate_sf)
  
  # normalize data of each cell based on computed size factors
  dge_norm <- t(t(dge) / size_factors) * mean(size_factors)
  
  # round to nearest number to generate normalized counts or scale using specified function
  if (norm_counts == TRUE) {
    round(dge_norm)
  }else if (!is.null(scale_fun)) {
    scale_fun(dge_norm)
  }
}

#' Test for differential expression
#' 
#' Test for differential expression across all perturbations. For each perturbation all positive
#' cells are selected and tested against cells carrying non-targeting negative controls.
#' 
#' @param perturb_status A data.frame containing binary perturbation information for each cell. Rows
#'   are perturbations (gRNAs or perturbed regulatory elements) and columns are cells, except for
#'   the first column, which has to contain the perturbation id.
#' @param gene_expr Matrix or data.frame containing digital gene expression data for all cells and
#'   genes to be tested. Rows are genes and columns are cells, with row names being gene ids and
#'   column names cell ids.
#' @param method Function to perform differential gene expression between perturbed and control
#'   cells.
#' @param parallel Triggers parallel computing using the \code{BiocParallel} package. This requires
#'   that a back-end was registered prior to executing the function.
test_differential_expression <- function(perturb_status, gene_expr, parallel = FALSE,
                                         method = c("MAST", "DEsingle", "LFC")) {
  
  # get specified DE method and attach required packages
  method <- match.arg(method)
  if (method %in% c("MAST", "DEsingle")) library(method, character.only = TRUE)
  
  # get lane barcodes of all cells
  cells <- colnames(gene_expr)
  i7_cells <- substr(cells, start = 1, stop = 8)
  cell_lanes <- data.frame(cell_barcode = cells, i7 = i7_cells, stringsAsFactors = FALSE)
  
  # get unique perturbation
  perts <- unique(perturb_status$VECTOR)
  names(perts) <- perts
  
  # perform differential gene expression test for each pertubation
  if (parallel == TRUE) {
    
    library(BiocParallel)
    output <- bplapply(perts, FUN = test_de, perturb_status = perturb_status, gene_expr = gene_expr,
                       cell_lanes = cell_lanes, method = method)
    
  }else{
    
    output <- lapply(perts, FUN = test_de, perturb_status = perturb_status, gene_expr = gene_expr,
                     cell_lanes = cell_lanes, method = method)
    
  }
  
  # convert output into one data.frame
  dplyr::bind_rows(output, .id = "perturbation")
  
}

#' DEsingle DE function
#' 
#' Perform differential expression analysis between perturbed and control cells using DEsingle.
#' 
#' @param seurat_obj Seurat object containing gene expression data and cell groupings.
de_DEsingle <- function(seurat_obj) {
  
  # normalize data using censored mean normalization
  counts <- GetAssayData(seurat_obj, slot = "raw.data")
  norm_counts <- normalize_cens_mean(counts, percentile = 0.9, norm_counts = TRUE)
  
  # get identity for each cell
  groups <- GetIdent(seurat_obj, uniq = FALSE)
  
  # detect DE genes
  de_results <- DEsingle(counts = norm_counts, group = groups)
  
  # reformat output
  data.frame(gene = rownames(de_results), de_results, stringsAsFactors = FALSE, row.names = NULL)
  
}

#' MAST DE function
#' 
#' Perform differential expression analysis between perturbed and control cells using MAST.
#' 
#' @param seurat_obj Seurat object containing gene expression data and cell groupings.
de_MAST <- function(seurat_obj) {
  
  # normalize data using censored mean normalization with log1p as scale function
  counts <- GetAssayData(seurat_obj, slot = "raw.data")
  norm_dat <- normalize_cens_mean(counts, percentile = 0.9, norm_counts = FALSE, scale_fun = log1p)
  seurat_obj <- SetAssayData(seurat_obj, assay.type = "RNA", slot = "data", new.data = norm_dat)

  # perform differential gene expression analysis using MAST
  de_results <- FindMarkers(seurat_obj, ident.1 = "pert", ident.2 = "ctrl", test.use = "MAST",
                            logfc.threshold = 0 , min.pct = 0)
  
  # reformat results
  data.frame(gene = rownames(de_results), de_results, stringsAsFactors = FALSE, row.names = NULL)
  
}

#' LFC DE function
#' 
#' Simple calculation of log-fold changes in gene expression between perturbed and control cells.
#' 
#' @param seurat_obj Seurat object containing gene expression data and cell groupings.
#' @param pseudocount Pseudocount to be added to transcript counts when calculating average gene
#'   expression.
de_LFC <- function(seurat_obj, pseudocount = 1) {
  
  # get raw counts and groups for each cell
  counts <- GetAssayData(seurat_obj, slot = "raw.data")
  groups <- GetIdent(seurat_obj, uniq = FALSE)
  
  # calculate average expression for each gene (with pseudocount)
  pert_avg <- rowMeans(counts[, groups == "pert"] + pseudocount)
  ctrl_avg <- rowMeans(counts[, groups == "ctrl"] + pseudocount)
  
  # merge into one data.frame and calculate lfc for each gene
  avg_genex <- data.frame(pert = pert_avg, ctrl = ctrl_avg)
  avg_genex$lfc <- log2(avg_genex$pert) - log2(avg_genex$ctrl)
  
  # reformat output
  data.frame(gene = rownames(avg_genex), avg_genex, stringsAsFactors = FALSE, row.names = NULL)
  
}
 
# helper functions ---------------------------------------------------------------------------------

# perform differential gene expression for one perturbation
test_de <- function(pert, perturb_status, gene_expr, cell_lanes, method, n_ctrl = 1000) {
  
  # get cells carrying perturbation pert
  pert_data <- perturb_status[perturb_status$VECTOR == pert, -1]
  cells_pert <- colnames(pert_data)[colSums(pert_data) > 0]
  
  # get 10x lane proportions for perturbed cells and calculate control cells per 10x lane
  i7_cells_pert <- cell_lanes[cell_lanes$cell_barcode %in% cells_pert, "i7"]
  i7_prop <- table(i7_cells_pert) / length(i7_cells_pert)
  i7_cells_ctrl <- round(i7_prop * n_ctrl)
  
  # randomly draw n (n_ctrl) control cells with same lane proportions as perturbed cells
  cell_lanes_ctrl <- cell_lanes[!cell_lanes$cell_barcode %in% cells_pert, ]  # not perturbed cells
  cells_ctrl <- lapply(names(i7_cells_ctrl), FUN = function(lane) {
    cells_lane <- cell_lanes_ctrl[cell_lanes_ctrl$i7 == lane, "cell_barcode"]
    sample(cells_lane, size = i7_cells_ctrl[lane], replace = FALSE)
  })
  cells_ctrl <- unlist(cells_ctrl)
  
  # create Seurat object containing gene expression for perturbed and sampled control cells
  seurat_obj <- CreateSeuratObject(raw.data = gene_expr[, c(cells_pert, cells_ctrl)],
                                   project = pert)
  
  # set identity (perturbed or control) for every cell
  groups <- c(rep("pert", length(cells_pert)), rep("ctrl", length(cells_ctrl)))
  seurat_obj <- SetIdent(seurat_obj, ident.use = groups)
  
  # get DE function for specified method
  de_function <- get(paste0("de_", method))
  
  # perform differential gene expression test. warnings and errors get reported and in case of
  # errors NULL is returned for that perturbation.
  tryCatch(
    withCallingHandlers({
      de_function(seurat_obj)
    }, warning = function(w) {
      message("For perturbation ", pert, ": ", w)
      invokeRestart("muffleWarning")
    }),
    error = function(e){
      message("For perturbation ", pert, ": ", e)
      return(NULL)
    })
}
