#!/usr/bin/env Rscript

# infer perturbation status of every cell based on simple vector transcript cutoff

# parse command line arguments ---------------------------------------------------------------------

library("optparse")

# create arguments list
option_list = list(
  make_option(c("-d", "--dge_file"), type = "character", default = NULL, 
              help = "Path to file containing digital expression (dge) data.",
              metavar = "character"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL,
              help = "Path to output file for perturbation status matrix", metavar = "character"),
  make_option(c("-v", "--vector_pattern"), type = "character", default = NULL,
              help = paste("Vector naming pattern (e.g. CROPseq_dCas9_DS).",
                            "Needed to extract vector expression data"),
              metavar = "character"),
  make_option(c("-t", "--min_txs"), type = "integer", default = NULL,
              help = "Minimum number (integer) of transcripts per vector to define perturbation",
              metavar = "integer")
)

# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser){
  
  if (is.null(opt[[arg]])){
    
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
    
  }
  
}

# check that all required parameters are provided
required_args <- c("dge_file", "outfile", "vector_pattern", "min_txs")
for (i in required_args) {
  
  check_required_args(i, opt = opt, opt_parser = opt_parser)
  
}

# define main function -----------------------------------------------------------------------------

# infer perturbation status of each cell based on simple transcript cutoff
infer_pert <- function(dge, vector_pattern, min_txs) {
  
  # filter for vector expression data
  vctrs <- grep(dge$GENE, pattern = vector_pattern)
  vctr_dge <- dge[vctrs, ]
  
  # make GENE to rownames and transpose
  rownames(vctr_dge) <- vctr_dge$GENE
  vctr_dge <- vctr_dge[, which(colnames(vctr_dge) != "GENE")]
  vctr_dge <- t(vctr_dge)
  
  # infer perturbation status of each cell
  pert <- (vctr_dge >= min_txs) * 1
  
  # add column for cell barcode
  pert <- cbind("cell_barcode" = rownames(pert), pert)
  rownames(pert) <- NULL
  
  return(pert)

}

# infer perturbation status of specified input file ------------------------------------------------

message("Loading DGE data...")
dge <- read.table(opt$dge_file, header = TRUE, stringsAsFactors = FALSE)

message("Inferring perturbation status...")
pert_status <- infer_pert(dge, vector_pattern = opt$vector_pattern, min_txs = opt$min_txs)

message("Writing to output file...")
write.table(pert_status, file = opt$outfile, quote = FALSE, sep = "\t")

message("Done!")
