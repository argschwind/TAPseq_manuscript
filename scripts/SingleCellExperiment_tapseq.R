## functions to generate SingleCellExperiment objects to store data from TAP-seq experiments

library(SingleCellExperiment)
library(methods)

#' Create a SingleCellExperiment object containing TAP-seq data
#' 
#' Both gene expression and perturbation data on each cell wil be stored as assays. Any additional
#' cell metadata can be added as well.
#' 
#' @param dge A data.frame or matrix containing digital gene expression data. Row names are assumed
#'   to be gene names and column names cell ids.
#' @param perturb_status A data.frame or matrix containing perturation status data. Row names are
#'   assumed to be perturbation ids and column names cell ids.
#' @param meta_data Optional DataFrame (or coercible to DataFrame) containing additional cell level
#'   meta data.
#' @param vector_pattern Optional character string providing pattern common to all vector ids. If
#'   specified this is used to remove any vector data from dge to only retain expression data on
#'   genes.
create_tapseq_sce <- function(dge, perturb_status, meta_data = NULL, vector_pattern = NULL) {
  
  # verify input
  if (any(colnames(dge) != colnames(perturb_status))) {
    stop("dge and perturb_status column names not identical! Not sorted correctly?", call. = FALSE)
  }
  
  # get rows in dge containing expression of genes and not vectors (if vector pattern is provided)
  if (!is.null(vector_pattern)) {
    gene_rows <- grep(rownames(dge), pattern = vector_pattern, invert = TRUE)
    dge <- dge[gene_rows, ]
  }
  
  # create sce object with gene expression
  if (is.null(meta_data)) {
    sce <- SingleCellExperiment(assays = list(counts = as.matrix(dge)))
  }else{
    sce <- SingleCellExperiment(assays = list(counts = as.matrix(dge)),
                                colData = DataFrame(meta_data))
  }
  
  # add perturbation status as internal metadata
  perturb_status <- DataFrame(t(as.matrix(perturb_status)))
  perturbations(sce) <- perturb_status
  
  return(sce)
  
}

# setters and getters for perturbation data --------------------------------------------------------

# define generics
setGeneric("perturbations", function(x) standardGeneric("perturbations"))
setGeneric("perturbations<-", function(x, value) standardGeneric("perturbations<-"))

# get perturbations from internal meta data
setMethod("perturbations", "SingleCellExperiment", function(x) {
  
  # get columns in internal meta data containing perturbations
  cd <- int_colData(x)
  pert_cols <- grepl(colnames(cd), pattern = "^PERT_.+$", perl = TRUE)
  perturb_status <- cd[, pert_cols]
  
  # remove internal identifier and return perturbation status
  colnames(perturb_status) <- sub("^PERT_", "", colnames(perturb_status), perl = TRUE)
  return(perturb_status)
  
})

# add perturbation status in internal meta data
setMethod("perturbations<-", c("SingleCellExperiment", "DataFrame"), function(x, value) {
  
  # verify input
  if (any(colnames(x) != rownames(value))) {
    stop("perturbation data rownames not the same as x column names. Not sorted correctly?",
         call. = FALSE)
  }
  
  # add prefix to perturbations, to separate from other internal colData
  colnames(value) <- paste0("PERT_", colnames(value))
  
  # get internal metadata and remove already existing perturbations
  ## FIND BETTER SOLUTION ASAP!!
  cd <- int_colData(x)
  pert_cols <- grepl(colnames(cd), pattern = "^PERT_.+$", perl = TRUE)
  cd <- cd[, !pert_cols]
  
  # add perturb status to internal meta data
  cd[, colnames(value)] <- value
  rownames(cd) <- rownames(value)
  int_colData(x) <- cd
  
  return(x)
  
})
