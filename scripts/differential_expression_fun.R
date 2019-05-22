## functions to perform differential gene expression analysis for enhancer TAP-seq ##

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
#' @param object A SingleCellExperiment object containing gene expression data and perturbation
#'   data as internal colData Any regular colData will be used as covariates if possible. Can be
#'   created by \code{create_tapseq_sce}.
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
#' @param p_adj_method Method to use for multiple testing correction. For more details see
#'   \code{\link[stats]{p.adjust}}.
test_differential_expression <- function(object, min_cells = 25, parallel = FALSE,
                                         method = c("MAST", "DEsingle", "LFC"),
                                         p_adj_method = c("fdr", "holm", "hochberg", "hommel",
                                                          "bonferroni", "BH", "BY", "none")) {
  
  # parse arguments and attach required packages (if needed)
  p_adj_method <- match.arg(p_adj_method)
  method <- match.arg(method)
  if (method %in% c("MAST", "DEsingle")) library(method, character.only = TRUE)
  
  # only retain perturbations with at least min_cells cells
  perts <- perturbations(object)
  cells_per_pert <- colSums(as.matrix(perts))
  perts <- perts[cells_per_pert >= min_cells]
  perturbations(object) <- perts
  
  # get lane barcodes of all cells
  cells <- colnames(object)
  i7_cells <- substr(cells, start = 1, stop = 8)
  cell_lanes <- data.frame(cell_barcode = cells, i7 = i7_cells, stringsAsFactors = FALSE)
  
  # get unique perturbation
  perts <- colnames(perturbations(object))
  names(perts) <- perts
  
  # perform differential gene expression test for each pertubation
  if (parallel == TRUE) {
    
    library(BiocParallel)
    output <- bplapply(perts, FUN = test_de, object = object, cell_lanes = cell_lanes,
                       method = method)
    
  }else{
    
    output <- lapply(perts, FUN = test_de, object = object, cell_lanes = cell_lanes,
                     method = method)
    
  }
  
  # convert output into one data.frame
  output <- dplyr::bind_rows(output, .id = "perturbation")
  
  # correct for multiple testing for all performed tests
  if (method != "LFC") {
    output$pval_adj_allTests <- p.adjust(output$pvalue, method = p_adj_method)
  }
  
  return(output)
  
}

#' DEsingle DE function
#' 
#' The first column in colData is expected to provide the perturbation status as a factor with two
#' levels, where the first level is non-perturbed cells. If this is not the case, the function will
#' attempt to convert the first column to a two level factor, with any risks that entails.
#' Any other columns in colData will ne ignored, because DEsingle currently does not support
#' covariates.
#'
#' @param pert_object SingleCellExperiment object containing gene expression data and cell groupings in
#'   colData. The perturbation to be tested is assumed to be the first column in colData!
de_DEsingle <- function(pert_object) {
  
  # normalize data using censored mean normalization
  counts <- assay(pert_object, "counts")
  norm_counts <- normalize_cens_mean(counts, percentile = 0.9, norm_counts = TRUE)
  
  # get identity for each cell and make sure it's a factor
  groups <- colData(pert_object)[, 1]
  groups <- as.factor(groups)
  levels(groups) <- c(0, 1)
  
  # detect DE genes
  de_results <- DEsingle(counts = norm_counts, group = groups)
  
  # reformat output
  data.frame(gene = rownames(de_results), de_results, stringsAsFactors = FALSE, row.names = NULL)
  
}

#' MAST DE function
#' 
#' Perform differential expression analysis between perturbed and control cells using MAST.
#'
#' The first column in colData is expected to provide the perturbation status as a factor with two
#' levels, where the first level is non-perturbed cells. If this is not the case, the function will
#' attempt to convert the first column to a two level factor, with any risks that entails.
#' Any other columns in colData will be used as covariates during modelling, but only the
#' perturbation will be tested for signficance.
#'
#' @param pert_object SingleCellExperiment object containing gene expression data and cell groupings
#'   in colData. The perturbation to be tested is assumed to be the first column in colData!
de_MAST <- function(pert_object) {
  
  # normalize data using censored mean normalization with log1p as scale function
  counts <- assay(pert_object, "counts")
  logcounts <- normalize_cens_mean(counts, percentile = 0.9, norm_counts = FALSE, scale_fun = log1p)
  assay(pert_object, "logcounts") <- logcounts
  
  # add some row- and colData expected by MAST (not used in DE tests)
  rowData(pert_object) <- cbind(rowData(pert_object), primerid = rownames(pert_object))
  colData(pert_object) <- cbind(colData(pert_object), wellKey = colnames(pert_object))
  
  # coerce pert_object to SingleCellAssay, since MAST requires that as input
  sca <- SceToSingleCellAssay(pert_object)
  
  # prepare colData for model fitting
  colnames(colData(sca))[1] <- "pert"  # make sure 1st column is called pert for coeff testing
  pert <- as.factor(colData(sca)$pert)
  levels(pert) <- c(0, 1)
  colData(sca)$pert <- pert
  
  # fit hurdle model
  vars <- colnames(colData(sca))
  vars_to_test <- vars[vars != "wellKey"]
  model <- paste0("~", paste(vars_to_test, collapse = "+"))
  zlm_fit <- zlm(as.formula(model), sca)
  
  # following command throws warning with MAST v1.8.1, so lrt test is performed separately
  # summary_zlm_fit <- summary(zlm_fit, logFC = TRUE, doLRT = "pert1")
  
  # get logFC change for the perturbation coefficient
  summary_zlm_fit <- summary(zlm_fit, logFC = TRUE)
  summary_dt <- summary_zlm_fit$datatable
  logFC <- summary_dt[contrast == "pert1" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)]
  
  # perform likelihood ratio test for the perturbation coefficient
  lrt <- as.data.frame(lrTest(zlm_fit, CoefficientHypothesis("pert1")))
  lrt_output <- data.frame(primerid = rownames(lrt), pvalue = lrt[, "hurdle.Pr(>Chisq)"],
                           row.names = rownames(lrt), stringsAsFactors = FALSE)
  
  # assemble and return output
  output <- merge(as.data.frame(logFC), lrt_output, by = "primerid")
  colnames(output) <- c("gene", "logFC", "ci_high", "ci_low", "pvalue")
  output[order(output$pvalue), ]
  
}

#' LFC DE function
#' 
#' Simple calculation of log-fold changes in gene expression between perturbed and control cells.
#' 
#' @param pert_object SingleCellExperiment object containing gene expression data (counts) and cell
#'   groupings in colData. The perturbation to be tested is assumed to be the first column in
#'   colData!
#' @param pseudocount Pseudocount to be added to transcript counts when calculating average gene
#'   expression.
de_LFC <- function(pert_object, pseudocount = 1) {
  
  # get raw counts and groups for each cell
  counts <- assay(pert_object, "counts")
  groups <- colData(pert_object)[,1]
  
  # calculate average expression for each gene (with pseudocount)
  pert_avg <- rowMeans(counts[, groups == 1] + pseudocount)
  ctrl_avg <- rowMeans(counts[, groups == 0] + pseudocount)
  
  # merge into one data.frame and calculate lfc for each gene
  avg_genex <- data.frame(pert = pert_avg, ctrl = ctrl_avg)
  avg_genex$lfc <- log2(avg_genex$pert) - log2(avg_genex$ctrl)
  
  # reformat output
  data.frame(gene = rownames(avg_genex), avg_genex, stringsAsFactors = FALSE, row.names = NULL)
  
}


# helper functions ---------------------------------------------------------------------------------

# perform differential gene expression for one perturbation
test_de <- function(pert, object, cell_lanes, method, n_ctrl = 1000) {
  
  # get cells carrying perturbation pert
  pert_data <- perturbations(object)[pert]
  cells_pert <- rownames(pert_data)[pert_data[, 1] > 0]
  
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
  
  # extract data for perturbed and sampled control cells
  pert_object <- object[, c(cells_pert, cells_ctrl)]
  
  # add perturbation status for "pert" as first column in colData, as this is needed for DE tests
  pert_status <- pert_data[colnames(pert_object), pert]
  cd <- cbind(pert = factor(pert_status), colData(pert_object), row.names = colnames(pert_object))
  levels(cd$pert) <- c(0, 1)
  cd <- droplevels(cd)  # remove any unused factor levels, as they might cause issues
  colData(pert_object) <- cd
  
  # get DE function for specified method
  de_function <- get(paste0("de_", method))
  
  # perform differential gene expression test. warnings and errors get reported and in case of
  # errors NULL is returned for that perturbation.
  tryCatch(
    withCallingHandlers({
      de_function(pert_object)
    }, warning = function(w) {
      message("For perturbation ", pert, ": ", w)
      invokeRestart("muffleWarning")
    }),
    error = function(e){
      message("For perturbation ", pert, ": ", e)
      return(NULL)
    })
  
}
