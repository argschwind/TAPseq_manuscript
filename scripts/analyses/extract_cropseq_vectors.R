#!/usr/bin/env Rscript

# script to extract CROP-seq vector sequences from aligned bam files. filters for reads that map to
# vector reference sequences in cells for which dge was extracted. cell barcodes for these cells are
# obtained from the dge summaryfile produced by Drop-seq tools based CROP-seq workflow.

# parse command line arguments =====================================================================

library(optparse)

# create arguments list
option_list = list(
  make_option(c("-b", "--bam_file"), type = "character", default = NULL,
              help = "Path to bam file containing aligned reads",
              metavar = "character"),
  make_option(c("-d", "--dge_summary_file"), type = "character", default = NULL,
              help = "Path to txt file containing dge summary", metavar = "character"),
  make_option(c("-p", "--vector_pattern"), type = "character", default = NULL,
              help = "Pattern for reads mapping to CROP-seq vectors", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "Output file", metavar = "character")
  
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
required_args <- c("bam_file", "dge_summary_file", "vector_pattern", "output_file")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

# define functions =================================================================================

# function to extract reads mapping to CROP-seq vectors and extract number of mismatches etc.
extract_vectors <- function(bam, cell_barcodes, output_file, vector_pattern) {
  
  # create cell barcode patterns and write to temporary file
  cbc_patterns <- paste0("XC:Z:", cell_barcodes)
  cbc_file <- tempfile(pattern = "cbc_patterns_", fileext = ".txt")
  write(cbc_patterns, file = cbc_file)
  
  # write heder to file
  header = c("cell_barcode", "umi_barcode", "aligned_vector_id", "mapping_quality", "mismatches",
             "mismatch_positions", "sequence")
  gz_outfile <- gzfile(output_file, open = "w")
  write(paste(header, collapse = "\t"), file = gz_outfile)
  close(gz_outfile)
  
  # bash command to filter bam file for CROP-seq vector reads from specified cells and to extract
  # mismatches, barcodes, and read sequence
  command <- paste("samtools view", bam, "| fgrep -f", cbc_file,
    "| fgrep", paste0("GE:Z:", vector_pattern),
    "| awk 'BEGIN { FS = \" \" }; {print $12\"\\t\"$20\"\\t\"$14\"\\t\"$5\"\\t\"$19\"\\t\"$13\"\\t\"$10}'",
    "| gzip >> ", output_file)
  
  # execute command
  system(command)
  
}

# extract CROP-seq vectors from provided bam file ==================================================

# read dge summary file and extract cell barcodes
dge_summary <- read.table(opt$dge_summary_file, header = TRUE, stringsAsFactors = FALSE)
cells <- dge_summary$cell_barcode

# extract CROP-seq vectors
extract_vectors(opt$bam_file, cell_barcodes = cells, output_file = opt$output_file,
                vector_pattern = opt$vector_pattern)
