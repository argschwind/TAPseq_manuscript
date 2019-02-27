#!/usr/bin/env Rscript

# infer perturbation status of every cell based on simple vector transcript cutoff

# parse command line arguments ---------------------------------------------------------------------

library("optparse")

# create arguments list
option_list = list(
  make_option(c("-i", "--infile"), type = "character", default = NULL,
              help = "Path to file containing digital expression (DGE) input data.",
              metavar = "character"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL,
              help = "Path to output file for perturbation status matrix", metavar = "character"),
  make_option(c("-v", "--vector_pattern"), type = "character", default = NULL,
              help = paste(
                     "Identifier common to all CROP-seq vectors, for instance a prefix or suffix",
                     "to the unique vector identifier. For instance in ID 'CROPseq_dCas9_DS_01'",
                     "the common identifier is 'CROPseq_dCas9_DS_' while '01' is the unique vector",
                     "identifier. This parameter is needed to extract vector expression data from",
                     "DGE data.", sep = " \n\t\t"
                     ),
              metavar = "character"),
  make_option(c("-m", "--min_txs"), type = "integer", default = NULL,
              help = "Minimum number (integer) of transcripts per vector to define perturbation",
              metavar = "integer"),
  make_option(c("-t", "--trim_id"), action = "store_true", default = FALSE,
              help = "Triggers removal of the vector pattern output vector IDs.")
)

# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  
  if (is.null(opt[[arg]])) {
    
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
    
  }
  
}

# check that all required parameters are provided
required_args <- c("infile", "outfile", "vector_pattern", "min_txs")
for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}

# define main function -----------------------------------------------------------------------------

# infer perturbation status of each cell based on simple transcript cutoff
infer_pert <- function(dge, vector_pattern, min_txs) {
  
  # filter for vector expression data
  vctrs <- grep(dge$GENE, pattern = vector_pattern)
  vctr_dge <- dge[vctrs, ]
  
  # move gene name to rownames
  rownames(vctr_dge) <- vctr_dge$GENE
  vctr_dge <- vctr_dge[, which(colnames(vctr_dge) != "GENE")]
  
  # infer perturbation status of each cell
  pert <- as.data.frame((vctr_dge >= min_txs) * 1)

  # add column for vector id
  pert <- cbind.data.frame("VECTOR" = rownames(pert), pert, stringsAsFactors = FALSE)
  rownames(pert) <- NULL
  
  return(pert)

}

# infer perturbation status of specified input file ------------------------------------------------

message("Loading DGE data...")
dge <- read.table(opt$infile, header = TRUE, stringsAsFactors = FALSE)

message("Inferring perturbation status...")
pert_status <- infer_pert(dge, vector_pattern = opt$vector_pattern, min_txs = opt$min_txs)

# remove vector pattern (common vector identifier) from output if specified
if (opt$trim_id == TRUE) {
  
  message("Trimming vector pattern from vector IDs...")
  pert_status$VECTOR <- sub(opt$vector_pattern, "", pert_status$VECTOR)
  
}

message("Writing to output file...")
write.table(pert_status, file = opt$outfile, quote = FALSE, sep = "\t", row.names = FALSE)

message("Done!")
