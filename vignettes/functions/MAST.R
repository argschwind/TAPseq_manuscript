library(MAST)

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

